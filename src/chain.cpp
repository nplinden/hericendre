#include "chain.h"
#include <Eigen/Sparse>
#include <cmath>
#include <deque>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

Chain::Chain(const std::string &path) : Chain(path.c_str()) {
}

Chain::Chain(const char *path) {
    // fmt::print("entering Chain\n");
    pugi::xml_document doc;
    doc.load_file(path);
    const pugi::xml_node chainxml = doc.child("depletion_chain");

    // INTIALIZATION STEP
    // std::cout << "Initialization:\n";
    for (pugi::xml_node nuclide = chainxml.child("nuclide"); nuclide;
         nuclide = nuclide.next_sibling("nuclide")) {
        nuclides_.push_back(std::make_shared<Nuclide>(Nuclide(nuclide)));

        nuclides_.back()->idInChain = nuclides_.size() - 1;
        for (pugi::xml_node decayNode = nuclide.child("decay"); decayNode;
             decayNode = decayNode.next_sibling("decay")) {
            decays_.push_back(
                std::make_shared<Decay>(Decay(decayNode, nuclides_.back())));

            std::string type = decays_.back()->type_;
            const double br = decays_.back()->branchingRatio_;
            if (Decay::SECONDARIES.find(type) != Decay::SECONDARIES.end()) {
                auto secVector = Decay::SECONDARIES.at(type);
                std::map<std::string, double> multiplicities;
                for (const auto &target: secVector) {
                    if (multiplicities.find(target) != multiplicities.end()) {
                        multiplicities[target] += br;
                    } else {
                        multiplicities[target] = br;
                    }
                }

                for (const auto &[k, v]: multiplicities) {
                    if (this->isIn(k)) {
                        decays_.push_back(std::make_shared<Decay>(
                            Decay(type, k, v, nuclides_.back())));
                    }
                }
            }
        }
        // for (pugi::xml_node reactionNode = nuclide.child("reaction");
        // reactionNode;
        //      reactionNode = reactionNode.next_sibling("reaction")) {
        //   if (reactionNode.attribute("type").value() != std::string("fission")) {
        //     reactions_.push_back(std::make_shared<NReaction>(
        //         NReaction(reactionNode, nuclides_.back())));
        //   }
        // }

        // for (pugi::xml_node nfyNode = nuclide.child("neutron_fission_yields");
        //      nfyNode; nfyNode = nfyNode.next_sibling("neutron_fission_yields")) {
        //   fissions_.push_back(
        //       std::make_shared<Fission>(Fission(nfyNode, nuclides_.back())));
        // }
    }
    // std::cout << "\t" << nuclides_.size() << " nuclides found\n";
    // std::cout << "\t" << decays_.size() << " decay reactions found\n";
    // std::cout << "\t" << reactions_.size() << " neutron reactions found\n";
    // std::cout << "\t" << fissions_.size() << " fissions found\n";

    // std::cout << "Binding decays\n";
    // BINDING STEP
    for (auto &dec: decays_) {
        dec->parent_->decays_.push_back(dec);
        if (!dec->targetName_.empty()) {
            dec->target_ = this->find(dec->targetName_);
            dec->target_->decaysUp_.push_back(dec);
        }
    }

    // std::cout << "Binding reaction\n";
    // for (auto &reac : reactions_) {
    //   reac->parent_->reactions_.push_back(reac);
    //   if (!reac->targetName_.empty())
    //     reac->target_ = this->find(reac->targetName_);
    //   reac->target_->reactionsUp_.push_back(reac);
    // }

    // for (auto& nuc : nuclides_){
    //   if (!nuc->nfyParent_.empty()){
    //     NuclidePtr target = this->find(nuc->nfyParent_) ;
    //     nuc->neutronFissionYields_ = target->neutronFissionYields_ ;
    //   }
    // }

    // std::cout << "Check...";
    // // COHERENCE CHECK
    // for (auto &n : nuclides_) {
    //   if (n->ndecay_ != n->decays_.size())
    //     throw std::invalid_argument(
    //         fmt::format("{} {} {}", n->name_, n->ndecay_, n->decays_.size()));
    //   if (n->nreac_ != n->reactions_.size())
    //     throw std::invalid_argument(
    //         fmt::format("{} {} {}", n->name_, n->nreac_,
    //         n->reactions_.size()));
    // }
    // std::cout << "done" << std::endl;
}

Chain::Chain() = default;

bool Chain::write(const char *path) {
    pugi::xml_document doc;
    auto root = doc.append_child("depletion_chain");

    for (const auto &nuclide: this->nuclides_) {
        nuclide->addNode(root);
    }
    return doc.save_file(path, "  ");
}

bool Chain::isIn(const std::string &name) const {
    return std::any_of(
        nuclides_.begin(),
        nuclides_.end(),
        [name](const NuclidePtr &p) { return name == p->name_; });
}

NuclidePtr Chain::find(int nucid) const {
    for (auto nuc: nuclides_) {
        if (nuc->zam_ == nucid) {
            return nuc;
        }
    }
    const std::string err_msg =
            fmt::format("Nuclide zam {} does not exist in the chain", nucid);
    throw std::invalid_argument(err_msg);
}

NuclidePtr Chain::find(const std::string &name) const {
    for (auto nuc: nuclides_) {
        if (nuc->name_ == name) {
            return nuc;
        }
    }
    std::string err_msg = fmt::format(
        "[find(std::string name)] Nuclide {} does not exist in the chain", name);
    throw std::invalid_argument(err_msg);
}

size_t Chain::nuclide_index(int nucid) const {
    for (size_t i = 0; i < nuclides_.size(); i++) {
        if (nuclides_[i]->zam_ == nucid)
            return i;
    }
    const std::string err_msg =
            fmt::format("Nuclide zam {} does not exist in the chain", nucid);
    throw std::invalid_argument(err_msg);
}

size_t Chain::nuclide_index(const std::string &name) const {
    for (size_t i = 0; i < nuclides_.size(); i++) {
        if (nuclides_[i]->name_ == name)
            return i;
    }
    std::string err_msg = fmt::format("[nuclide_index(std::string name)] Nuclide "
                                      "{} does not exist in the chain",
                                      name);
    throw std::invalid_argument(err_msg);
}

Eigen::SparseMatrix<double> Chain::decayMatrix() const {
    std::vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(10000);
    for (size_t inuc = 0; inuc < this->nuclides_.size(); inuc++) {
        const NuclidePtr nuc = this->nuclides_[inuc];
        triplets.emplace_back(inuc, inuc, -nuc->dconst_);
        for (const auto &d: nuc->decays_) {
            if (!d->targetName_.empty()) {
                size_t jnuc = this->nuclide_index(d->targetName_);
                triplets.emplace_back(jnuc, inuc, nuc->dconst_ * d->branchingRatio_);
            }
        }
    }
    Eigen::SparseMatrix<double> M(this->nuclides_.size(), this->nuclides_.size());
    M.setFromTriplets(triplets.begin(), triplets.end());
    M.makeCompressed();
    return M;
}

void Chain::dfs(const size_t &nucid, std::vector<bool> &visited) {
    if (visited[nucid])
        return;
    visited[nucid] = true;

    std::vector<DecayPtr> decays = nuclides_[nucid]->decays_;
    std::vector<size_t> neighbours;

    for (const DecayPtr &decay: nuclides_[nucid]->decays_) {
        neighbours.push_back(decay->target_->idInChain);
    }
    for (const auto neighboursId: neighbours) {
        this->dfs(neighboursId, visited);
    }
}

std::vector<std::string> Chain::reachable(const std::string &nucname) {
    const size_t n = nuclides_.size();
    std::vector<bool> visited;
    for (size_t i = 0; i < n; i++)
        visited.push_back(false);

    const size_t initialId = this->find(nucname)->idInChain;
    this->dfs(initialId, visited);

    std::vector<std::string> names;
    for (size_t i = 0; i < n; i++) {
        if (visited[i])
            names.push_back(nuclides_[i]->name_);
    }
    return names;
}

bool Chain::topological_sort() {
    std::vector<NuclidePtr> sorted;
    // Building the vector of incoming degrees
    std::vector<size_t> incoming_degrees;
    incoming_degrees.reserve(this->nuclides_.size());
    for (const auto &nuclide: this->nuclides_) {
        incoming_degrees.push_back(nuclide->decaysUp_.size());
        // fmt::print("{:8} {}\n", nuclide->name_, nuclide->decaysUp_.size());
    }

    // Initialize the queue with orphan nuclides
    std::deque<NuclidePtr> queue;
    for (size_t inuc = 0; inuc < incoming_degrees.size(); inuc++) {
        if (incoming_degrees[inuc] == 0) {
            queue.push_back(this->nuclides_[inuc]);
        }
    }

    while (!queue.empty()) {
        NuclidePtr front = queue.front();
        sorted.push_back(front);
        for (const auto &decay: front->decays_) {
            if (decay->hasTarget_) {
                const size_t targetId = decay->target_->idInChain;
                incoming_degrees[targetId]--;
                if (incoming_degrees[targetId] == 0)
                    queue.push_back(this->nuclides_[targetId]);
            }
        }
        queue.pop_front();
    }

    if (sorted.size() != this->nuclides_.size()) {
        return false;
    } else {
        this->nuclides_ = sorted;
        for (size_t i = 0; i < this->nuclides_.size(); i++)
            this->nuclides_[i]->idInChain = i;

        return true;
    }
}

void Chain::tweak_dconst() {
    std::vector<double> dconsts;
    std::map<double, std::vector<NuclidePtr> > duplicates;
    for (const auto &nuclide: this->nuclides_) {
        duplicates[nuclide->dconst_].push_back(nuclide);
    }
    for (auto const &[key, val]: duplicates) {
        if (val.size() > 1) {
            for (size_t i = 0; i < val.size(); i++) {
                val[i]->dconst_ *= pow((1 + 1e-14), i);
            }
        }
    }
}

void Chain::dump_matrix(const std::string &path) const {
    Eigen::MatrixXd dMat(this->decayMatrix());
    // fmt::print("size = [{}, {}]\n", dMat.rows(), dMat.cols());

    auto out = fmt::output_file(path);
    out.print(",");
    for (size_t inuc = 0; inuc < this->nuclides_.size(); inuc++) {
        out.print("{}", this->nuclides_[inuc]->name_);
        if (inuc < this->nuclides_.size() - 1)
            out.print(",");
    }
    out.print("\n");
    for (long int row = 0; row < dMat.rows(); row++) {
        out.print("{},", this->nuclides_[row]->name_);
        for (long int col = 0; col < dMat.cols(); col++) {
            out.print("{}", dMat(row, col));
            if (col < dMat.cols() - 1)
                out.print(",");
        }
        out.print("\n");
    }
    out.close();
}

std::vector<std::string> Chain::name_vector() const {
    std::vector<std::string> vec ;
    vec.reserve(this->nuclides_.size());
    for (const auto& nuclide: this->nuclides_)
        vec.push_back(nuclide->name_);
    return vec;
}

std::vector<double> Chain::dconst_vector() const {
    std::vector<double> vec;
    vec.reserve(this->nuclides_.size());
    for (const auto& nuclide: this->nuclides_)
        vec.push_back(nuclide->dconst_);
    return vec;
}

// void Chain::removeNuclide(std::string nuc) {}

// void Chain::restrict(std::vector<std::string> allowed) {
// std::vector<NuclidePtr> newNuclides ;
// for (auto& nuc : nuclides_){
//   if (std::find(allowed.begin(), allowed.end(), nuc->name_) !=
//   allowed.end()){
//     newNuclides.push_back(nuc) ;
//   }
// }

// std::vector<DecayPtr> newDecays ;
// for (auto& decay : decays_){
//   bool parentOk = std::find(allowed.begin(), allowed.end(),
//   decay->parentName_) != allowed.end() ; bool targetOk =
//   std::find(allowed.begin(), allowed.end(), decay->targetName_) !=
//   allowed.end() ; if (parentOk && targetOk){
//     newDecays.push_back(decay) ;
//   }
// }
// }
