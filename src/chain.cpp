#include "chain.h"
#include <Eigen/Sparse>
#include <cmath>
#include <deque>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <highfive/H5Easy.hpp>

Chain::Chain(const std::string &path) : Chain(path.c_str())
{
}

Chain::Chain(const char *path)
{
    pugi::xml_document doc;
    pugi::xml_parse_result results = doc.load_file(path);
    if (!results)
    {
        std::string err_msg = fmt::format("Failed to load depletion chain file {} ({})", path, results.description());
        throw std::runtime_error(err_msg);
    }

    const pugi::xml_node root = doc.child("depletion_chain");
    if (!root)
    {
        std::string err_msg = fmt::format("File {} is not a valid depletion chain file", path);
        throw std::runtime_error(err_msg);
    }

    /*
    The chain is built in two steps to allow for nuclide to be referenced
    before they are defined in the xml file.
    */

    // Initialization step
    // This step adds all nuclides and decay reactions to the chain object.
    size_t nuclide_count = std::distance(root.children("nuclide").begin(), root.children("nuclide").end());
    nuclides_.reserve(nuclide_count);
    for (pugi::xml_node nuclide : root.children("nuclide"))
    {

        nuclides_.push_back(std::make_shared<Nuclide>(Nuclide(nuclide)));
        nuclides_.back()->idInChain = nuclides_.size() - 1;

        for (pugi::xml_node decayNode : nuclide.children("decay"))
        {
            decays_.push_back(std::make_shared<Decay>(Decay(decayNode, nuclides_.back())));

            if (decays_.back()->hasSecondaries())
            {
                decays_.push_back(std::make_shared<Decay>(decays_.back()->getSecondaries()));
            }
        }
    }

    // Binding step
    /*
    The binding step connects decay targets to the corresponding decay object.
    */
    for (auto &dec : decays_)
    {
        dec->parent_->decays_.push_back(dec);
        if (!dec->targetName_.empty())
        {
            dec->target_ = this->find(dec->targetName_);
            dec->target_->decaysUp_.push_back(dec);
        }
    }
}

Chain::Chain() = default;

bool Chain::write(const char *path)
{
    pugi::xml_document doc;
    auto root = doc.append_child("depletion_chain");

    for (const auto &nuclide : this->nuclides_)
    {
        nuclide->addNode(root);
    }
    return doc.save_file(path, "  ");
}

bool Chain::contains(const std::string &name) const
{
    for (const auto &p : nuclides_)
    {
        if (name == p->name_)
        {
            return true;
        }
    }
    return false;
}

NuclidePtr Chain::find(int nucid) const
{
    for (auto nuc : nuclides_)
    {
        if (nuc->zam_ == nucid)
        {
            return nuc;
        }
    }
    const std::string err_msg =
        fmt::format("Nuclide zam {} does not exist in the chain", nucid);
    throw std::invalid_argument(err_msg);
}

NuclidePtr Chain::find(const std::string &name) const
{
    for (auto nuc : nuclides_)
    {
        if (nuc->name_ == name)
        {
            return nuc;
        }
    }
    std::string err_msg = fmt::format(
        "[find(std::string name)] Nuclide {} does not exist in the chain", name);
    throw std::invalid_argument(err_msg);
}

size_t Chain::nuclide_index(int nucid) const
{
    for (size_t i = 0; i < nuclides_.size(); i++)
    {
        if (nuclides_[i]->zam_ == nucid)
            return i;
    }
    const std::string err_msg =
        fmt::format("Nuclide zam {} does not exist in the chain", nucid);
    throw std::invalid_argument(err_msg);
}

size_t Chain::nuclide_index(const std::string &name) const
{
    for (size_t i = 0; i < nuclides_.size(); i++)
    {
        if (nuclides_[i]->name_ == name)
            return i;
    }
    std::string err_msg = fmt::format("[nuclide_index(std::string name)] Nuclide "
                                      "{} does not exist in the chain",
                                      name);
    throw std::invalid_argument(err_msg);
}

Eigen::SparseMatrix<double> Chain::decayMatrix() const
{

    // eigen triplets are (row, col, value)
    const size_t n = nuclides_.size();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n + 3 * n); // rough estimate

    for (size_t inuc = 0; inuc < n; inuc++)
    {
        const NuclidePtr &nuc = nuclides_[inuc];
        triplets.emplace_back(inuc, inuc, -nuc->dconst_);

        for (const auto &d : nuc->decays_)
        {
            if (d->target_)
            {
                const size_t jnuc = d->target_->idInChain;
                triplets.emplace_back(jnuc, inuc, nuc->dconst_ * d->branchingRatio_);
            }
        }
    }
    Eigen::SparseMatrix<double> M(n, n);
    M.setFromTriplets(triplets.begin(), triplets.end());
    M.makeCompressed();
    return M;
}

void Chain::dfs(const size_t &nucid, std::vector<bool> &visited)
{
    if (visited[nucid])
        return;
    visited[nucid] = true;

    std::vector<DecayPtr> decays = nuclides_[nucid]->decays_;
    std::vector<size_t> neighbours;

    for (const DecayPtr &decay : nuclides_[nucid]->decays_)
    {
        neighbours.push_back(decay->target_->idInChain);
    }
    for (const auto neighboursId : neighbours)
    {
        this->dfs(neighboursId, visited);
    }
}

std::vector<std::string> Chain::reachable(const std::string &nucname)
{
    const size_t n = nuclides_.size();
    std::vector<bool> visited;
    for (size_t i = 0; i < n; i++)
        visited.push_back(false);

    const size_t initialId = this->find(nucname)->idInChain;
    this->dfs(initialId, visited);

    std::vector<std::string> names;
    for (size_t i = 0; i < n; i++)
    {
        if (visited[i])
            names.push_back(nuclides_[i]->name_);
    }
    return names;
}

bool Chain::topological_sort()
{
    std::vector<NuclidePtr> sorted;
    // Building the vector of incoming degrees
    std::vector<size_t> incoming_degrees;
    incoming_degrees.reserve(this->nuclides_.size());
    for (const auto &nuclide : this->nuclides_)
    {
        incoming_degrees.push_back(nuclide->decaysUp_.size());
        // fmt::print("{:8} {}\n", nuclide->name_, nuclide->decaysUp_.size());
    }

    // Initialize the queue with orphan nuclides
    std::deque<NuclidePtr> queue;
    for (size_t inuc = 0; inuc < incoming_degrees.size(); inuc++)
    {
        if (incoming_degrees[inuc] == 0)
        {
            queue.push_back(this->nuclides_[inuc]);
        }
    }

    while (!queue.empty())
    {
        NuclidePtr front = queue.front();
        sorted.push_back(front);
        for (const auto &decay : front->decays_)
        {
            if (decay->hasTarget_)
            {
                const size_t targetId = decay->target_->idInChain;
                incoming_degrees[targetId]--;
                if (incoming_degrees[targetId] == 0)
                    queue.push_back(this->nuclides_[targetId]);
            }
        }
        queue.pop_front();
    }

    if (sorted.size() != this->nuclides_.size())
    {
        return false;
    }
    else
    {
        this->nuclides_ = sorted;
        for (size_t i = 0; i < this->nuclides_.size(); i++)
            this->nuclides_[i]->idInChain = i;

        return true;
    }
}

void Chain::tweak_dconst()
{
    // Group nuclides by decay constant
    std::map<double, std::vector<NuclidePtr>> duplicates;
    for (const auto &nuclide : nuclides_)
    {
        duplicates[nuclide->dconst_].push_back(nuclide);
    }

    // Add small perturbations to decay constants
    for (auto const &[_, nuclides] : duplicates)
    {
        if (nuclides.size() > 1)
        {
            for (size_t i = 0; i < nuclides.size(); i++)
            {
                nuclides[i]->dconst_ *= std::pow((1 + 1e-7), i);
            }
        }
    }
}

void Chain::save(const std::string &path) const
{
    Eigen::SparseMatrix<double> dMat(this->decayMatrix());

    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<size_t> rows;
    std::vector<size_t> cols;
    std::vector<double> values;

    triplets.reserve(dMat.nonZeros());

    for (int k = 0; k < dMat.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(dMat, k); it; ++it)
        {
            rows.push_back(it.row());
            cols.push_back(it.col());
            values.push_back(it.value());
        }
    }
    H5Easy::File file(path, H5Easy::File::Overwrite);
    H5Easy::dump(file, "ROW", rows);
    H5Easy::dump(file, "COL", cols);
    H5Easy::dump(file, "VALUES", values);
}

std::vector<std::string> Chain::name_vector() const
{
    std::vector<std::string> vec;
    vec.reserve(this->nuclides_.size());
    for (const auto &nuclide : this->nuclides_)
        vec.push_back(nuclide->name_);
    return vec;
}

std::vector<double> Chain::dconst_vector() const
{
    std::vector<double> vec;
    vec.reserve(this->nuclides_.size());
    for (const auto &nuclide : this->nuclides_)
        vec.push_back(nuclide->dconst_);
    return vec;
}
