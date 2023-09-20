from __future__ import annotations
import numpy as np
import xml.etree.ElementTree as ET
import scipy as sp
import networkx as nx
from itertools import groupby
from hericendre import Nuclide


class Chain(nx.MultiDiGraph):
    def __init__(self, initiator=None):
        super(Chain, self).__init__()
        if isinstance(initiator, nx.MultiDiGraph):
            super(Chain, self).__init__(initiator)

    @classmethod
    def from_xml(cls, path: str) -> Chain:
        chain = cls()
        f = ET.parse(path)
        sf_nfy_parents = {}
        sf_nfys = {}
        sf_nfy_br = {}
        nfy_parents = {}
        nfy_Q = {}
        nfys = {}
        for isotnode in f.findall("nuclide"):
            parent = isotnode.get("name")
            hl = isotnode.get("half_life")
            dconst = np.log(2) / float(hl) if hl is not None else 0.
            decay_energy = isotnode.get("decay_energy", 0.)

            if parent not in chain:
                chain.add_node(parent, dconst=dconst, decay_energy=decay_energy, zae=Nuclide(parent).nucid)
            else:
                chain.nodes[parent]["dconst"] = dconst
                chain.nodes[parent]["decay_energy"] = decay_energy
                chain.nodes[parent]["zae"] = Nuclide(parent).nucid

            for reaction in isotnode.findall("reaction"):
                if reaction.get("type") == "fission":
                    nfy = isotnode.find("neutron_fission_yields")
                    Q = reaction.get("Q", default=0.)
                    if (nfy_parent := nfy.get("parent")) is not None:
                        nfy_parents[parent] = nfy_parent
                        nfy_Q[parent] = Q
                    else:
                        nfys[parent] = {}
                        for fy in nfy.findall("fission_yields"):
                            products = fy.find("products").text.split()
                            data = [float(val) for val in fy.find("data").text.split()]
                            energy = fy.get("energy")
                            nfys[parent][energy] = {
                                "products": products,
                                "data": data
                            }
                            for pf, br in zip(products, data):
                                chain.add_edge(
                                    parent,
                                    pf,
                                    key=energy,
                                    event="reaction",
                                    type="Fission",
                                    br=br,
                                    Q=Q
                                )
                else:
                    chain.add_edge(
                        parent,
                        reaction.get("target", default="void"),
                        event="reaction",
                        type=reaction.get("type"),
                        br=float(reaction.get("branching_ratio", default=1.)),
                        Q=reaction.get("Q", default=0.)
                    )

            for decaynode in isotnode.findall("decay"):
                if decaynode.get("type") == "sf":
                    decay_br = float(decaynode.get("branching_ratio", default=1.))
                    nfy = isotnode.find("neutron_fission_yields")
                    try:
                        if (nfy_parent := nfy.get("parent")) is not None:
                            sf_nfy_parents[parent] = nfy_parent
                            sf_nfy_br[parent] = decay_br
                        else:
                            sf_nfys[parent] = {}
                            fy = sorted([e for e in nfy.findall("fission_yields")],
                                        key=lambda x: float(x.get("energy")))[0]
                            products = fy.find("products").text.split()
                            data = [float(val) for val in fy.find("data").text.split()]
                            sf_nfys[parent] = {
                                "products": products,
                                "data": data
                            }
                            for pf, br in zip(products, data):
                                chain.add_edge(
                                    parent,
                                    pf,
                                    key="sf",
                                    event="decay",
                                    type="Fission",
                                    br=br * decay_br,
                                )
                    except AttributeError:
                        chain.add_edge(
                            parent,
                            "void",
                            key="sf",
                            event="decay",
                            type="Fission",
                            br=decay_br,
                        )

                else:
                    chain.add_edge(
                        parent,
                        decaynode.get("target", default="void"),
                        event="decay",
                        type=decaynode.get("type"),
                        br=float(decaynode.get("branching_ratio", default=1.))
                    )

        for parent, substitute in nfy_parents.items():
            for energy in nfys[substitute]:
                products = nfys[substitute][energy]["products"]
                data = nfys[substitute][energy]["data"]
                for pf, br in zip(products, data):
                    chain.add_edge(
                        parent,
                        pf,
                        key=energy,
                        event="reaction",
                        type="Fission",
                        br=br,
                        Q=nfy_Q[parent]
                    )
        for parent, substitute in sf_nfy_parents.items():
            products = sf_nfys[substitute]["products"]
            data = sf_nfys[substitute]["data"]
            for pf, br in zip(products, data):
                chain.add_edge(
                    parent,
                    pf,
                    key=0,
                    event="decay",
                    type="Fission",
                    br=br * sf_nfy_br[parent],
                )

        chain.nodes["void"]["dconst"] = 0.
        chain.nodes["void"]["decay_energy"] = 0.
        chain.nodes["void"]["zae"] = -1
        return chain

    @classmethod
    def from_mdl(cls, path: str) -> Chain:
        chain = cls()
        branchelines = []
        fissionlines = []
        with open(path, "r") as f:
            fissionflag = False
            for line in f.readlines():
                if cls.iscomment(line):
                    continue
                else:
                    if len(line.split()) == 1:
                        fissionflag = True
                    if not fissionflag:
                        branchelines.append(line.split())
                    else:
                        fissionlines.append(line.split())

        daughter = ""
        for line in branchelines:
            if line[2].isdigit():  # A daughter node
                daughter = line[0]
                if daughter not in chain:
                    chain.add_node(daughter)
                attrs = {daughter: {"dconst": float(line[1]), }}
                nx.set_node_attributes(chain, attrs)
            else:
                parent = line[0]
                reaction = line[1]
                rtype = "decay" if "DRTYP" in reaction else "reaction"
                br = float(line[2])
                chain.add_edge(parent, daughter, event=rtype, name=reaction, br=br)
        return chain

    def treat_equal_dconst(self, epsilon: float = 1e-12):
        c = {n: self.nodes[n]["dconst"] for n in self.nodes if self.nodes[n]["dconst"] != 0}

        duplicates = {}
        for k, v in c.items():
            duplicates[v] = duplicates.get(v, []) + [k]
        duplicates = {k: v for k, v in duplicates.items() if len(v) > 1}

        lens = {k: len(v) for k, v in duplicates.items()}
        for dconst, isotopes in duplicates.items():
            for i, isot in enumerate(isotopes):
                self.nodes[isotopes[i]]["dconst"] *= (1 + epsilon)**i


    def ascendance(self, isot):
        view = nx.graphviews.reverse_view(self)
        return {daughter: dict(view[isot][daughter]) for daughter in view[isot]}

    def descendance(self, isot):
        return {daughter: dict(self[isot][daughter]) for daughter in self[isot]}

    def dconst(self, isot):
        return self.nodes[isot]["dconst"]

    def stable(self, isot):
        return self.dconst(isot) == 0.

    def halflife(self, isot):
        try:
            return 1 / self.nodes[isot]["dconst"]
        except ZeroDivisionError:
            raise ZeroDivisionError(f"Isotope {isot} is stable.")

    @property
    def acyclical(self):
        return nx.is_directed_acyclic_graph(self)

    def zero_power(self, inplace=True):
        if not inplace:
            reduced = Chain(self)
            reduced.remove_edges_from(
                [edge for edge in reduced.edges if reduced.edges[edge]["event"] != "decay"]
            )
            return reduced
        else:
            self = Chain(self)
            self.remove_edges_from(
                [edge for edge in self.edges if self.edges[edge]["event"] != "decay"]
            )

    def remove_void(self, inplace=True):
        if not inplace:
            reduced = Chain(self)
            reduced.remove_edges_from(
                [edge for edge in reduced.edges if "void" in edge]
            )
            reduced.remove_node("void")
            return reduced
        else:
            self.remove_edges_from(
                [edge for edge in self.edges if "void" in edge]
            )
            self.remove_node("void")

    def topological_order(self):
        return list(nx.topological_sort(self))

    def build_decay_matrix(self, order: str = "topo"):
        reduced = self.zero_power(inplace=False) #.remove_void(inplace=False)
        if "Es258" in reduced: # A hack for removing the Es258 -> Es258 loop of ENDF-B/VII.1 chains
            reduced.remove_node("Es258")
        assert reduced.acyclical

        if order == "file":
            isotopes = list(reduced.nodes)
        elif order == "alpha":
            isotopes = sorted(list(reduced.nodes))
        elif order == "zae":
            isotopes = sorted(
                list(reduced.nodes), key=lambda x: reduced.nodes[x]["zae"]
            )
        elif order == "topo":
            isotopes = reduced.topological_order()
        else:
            raise ValueError("Unknow order type.")
        # matrix = np.zeros((len(isotopes), len(isotopes)))
        nisot = len(isotopes)
        matrix = sp.sparse.dok_array((nisot, nisot))
        for i, parent in enumerate(isotopes):
            matrix[i, i] = -reduced.dconst(parent)
            for daughter in reduced[parent]:
                j = isotopes.index(daughter)
                for _, edge in reduced[parent][daughter].items():
                    matrix[j, i] = reduced.dconst(parent) * edge["br"]
        return isotopes, matrix

    @staticmethod
    def iscomment(line_):
        comment_symb = ("//", "/*", "#")
        b = line_.isspace() or line_.strip().startswith(comment_symb)
        return b
