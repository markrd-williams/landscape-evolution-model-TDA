import numpy as np
import gudhi.wasserstein
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

"""
This file contains the matching approximate vineyard calculations.
This also contains the code for performing statistics on the approximate vineyards.
For the real vineyards there (will be) another file containing the code for that.
"""


class DiscreteVineyard:
    def __init__(self, dgms):
        """ Form a discrete vineyard from the diagram time series input.

        arguments:
        - dgms: A time series of the H0 component or augmented H0 component of persistence diagrams

        computes graph in the form of an adjacency list.
        Vertices are of the form (t,i) where dgms[t,i] is the homology class it represents
        Edges of the form (t, i), (t+1, j) if (i,j) is an edge in a Wasserstein matching between dgms[t], dgms[t+1]
        """
        if len(dgms[0][0]) > 2:
            # In this case the diagrams must be augmented
            self.augmented = True
        else:
            self.augmented = False

        self.dgms = dgms
        self._graph()
        self._vines()

    def _graph(self):
        graph = {}
        # i is the timestep the prev_diagram is from
        for t, (prev_dgm, cur_dgm) in enumerate(zip(self.dgms, self.dgms[1:])):
            _, matching = gudhi.wasserstein.wasserstein_distance(prev_dgm,
                                                                 cur_dgm,
                                                                 matching=True,
                                                                 internal_p=2,
                                                                 order=2)
            # i and j are indices of the components matched
            for [i,j] in matching:
                if i == -1:
                    # In this case we are matching something from cur_dgm to the diagonal of prev_dgm
                    # This means that we don't care about it in this stage
                    continue
                elif j != -1:
                    graph[(t, i)] = [(t+1, j)]
                    graph[(t+1, j)] = [] # This is to ensure there always exists a node for each one we reference
                else:
                    graph[(t, i)] = []

        self.graph = graph

    def _vines(self):
        """Output a list of the vines (i.e. connected components) of the discrete vineyard"""
        vines = []
        seen = {}
        for key in self.graph:
            seen[key] = False
            # sorting by layer, this way we always find the first node on a vine
        for root in sorted(self.graph.keys()):
            if not seen[root]:
                # Do a DFS on the unseen node to get vine
                vine = [root]
                stack = [root]
                while stack:
                    v = stack.pop()
                    if not seen[v]:
                        seen[v] = True
                        for w in self.graph[v]:
                            stack.append(w)
                            vine.append(w)
                vines.append(vine)
        self.vines = vines

    def stats(self):
        """Output two statistics, which can be used to compute other statistics.

        Returns:
        - born_at_layer: List of lists. Each entry is for a time timestep. Each entry is a list of persistences of new components born at that time
        - death_at_layer: Same but for death
        """
        born_at_layer = [[] for _ in range(len(self.dgms))]
        died_at_layer = [[] for _ in range(len(self.dgms))]
        for vine in self.vines:
            born_t, born_i = vine[0]
            death_t, death_i = vine[-1]

            born_b, born_d, *_ = self.dgms[born_t][born_i]
            death_b, death_d, *_ = self.dgms[death_t][death_i]

            born_at_layer[born_t].append( abs(born_b - born_d) )
            died_at_layer[death_t].append( abs(death_b - death_d ) )
        return born_at_layer, died_at_layer

    # TODO: Would this be faster from vines?
    #      That way can add longer segments.
    def plot(self, ax=None, xlim=[-np.inf, np.inf], ylim=[-np.inf, np.inf], zlim=[-np.inf, np.inf]):
        if not ax:
            ax = plt.axes(projection="3d")

        # Since line collection does not keep track of the extent of the image we want to display,
        # we keep track of bounds during computation
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf
        zmin = np.inf
        zmax = -np.inf

        segments = []
        colours = []
        for (layer, comp_index), nodes in self.graph.items():
            if layer > zlim[1]:
                break
            if layer < zlim[0]:
                continue
            for layer2, comp2_index in nodes:
                if layer2 > zlim[1]:
                    continue
                c1_birth, c1_death, *_ = self.dgms[layer][comp_index]
                c2_birth, c2_death, *_ = self.dgms[layer2][comp2_index]
                if (c1_birth < xlim[0] or c1_birth > xlim[1] or
                    c1_death < ylim[0] or c1_death > ylim[1] or
                    c2_birth < xlim[0] or c2_birth > xlim[1] or
                    c2_death < ylim[0] or c2_birth > ylim[1]):
                    continue
                segments.append([[c1_birth, c1_death, layer], [c2_birth, c2_death, layer2]])
                colours.append(cmap(255 * layer // len(dgms) ))

                xmin = min(xmin, c1_birth, c2_birth)
                xmax = max(xmax, c1_birth, c2_birth)
                ymin = min(ymin, c1_death, c2_death)
                if not np.isinf(c1_death):
                    ymax = max(ymax, c1_death)
                if not np.isinf(c2_death):
                    ymax = max(ymax, c2_death)
                zmin = min(zmin, layer, layer2)
                zmax = max(zmax, layer, layer2)

        lc = Line3DCollection(segments, colors=colours)
        ax.add_collection3d(lc)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_zlim(zmin, zmax)
        return ax

    def average_vine_length(self):
        max_time = len(self.dgms)
        values = np.zeros(max_time)
        number = np.zeros(max_time)
        for v in self.vines:
            born_time = v[0][0]
            values[born_time] += len(v)
            number[born_time] += 1
        return values / number

    def average_peak_travel(self):
        if not self.augmented:
            raise Exception("Diagrams input not augmented, cannot find average peak travel")
        avg = 0
        for vine in self.vines:
            dist = 0
            for (t,i),(tp,ip) in zip(vine, vine[1:]):
                disp = self.dgms[t][i][2:] - self.dgms[tp][ip][2:]
                dist += math.sqrt(disp[0] ** 2 + disp[1] ** 2)
            avg += dist
        return avg / len(vines)
