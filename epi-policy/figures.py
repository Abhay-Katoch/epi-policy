import matplotlib.pyplot as plt
import networkx as nx

############################################################
# Figure 1
############################################################

# Use default matplotlib mathtext renderer instead of LaTeX
plt.rcParams['text.usetex'] = False


nodes = ["S", "E", "P", "I", "A", "R", "D"]

# Instantiate nodes for DiGraph
G = nx.DiGraph()
G.add_nodes_from(nodes)

# Add edges with LaTeX formatted labels for each node. 
edges = [
    ('S', 'E', {'label': r'$\lambda_t$'}),
    ('E', 'P', {'label': r'$\sigma$'}),
    ('P', 'I', {'label': r'$\delta(1-\rho)$'}),
    ('P', 'A', {'label': r'$\delta \cdot \rho$'}),
    ('I', 'R', {'label': r'$\gamma_I$'}),
    ('A', 'R', {'label': r'$\gamma_A$'}),
    ('I', 'D', {'label': r'$\gamma \cdot r$'}),
]

# Manually set the positions for each node
positions = {
    'S': (0, 0.5),
    'E': (1, 0.5),
    'P': (2, 0.5),
    'I': (3, 1),
    'A': (3, 0),
    'R': (4, 0),
    'D': (4, 1)
}

# Add the edges between each node for the directed graph
G.add_edges_from((u, v, attr) for (u, v, attr) in edges)

# Draw the graph
nx.draw(
    G,
    positions,
    with_labels = True,
    node_size=2000, 
    node_color = "lightgrey", 
    font_size=10, 
    font_weight = "bold", 
    edge_color = "black"
)

# Draw edge labels from the list of edges and directed graph G.
nx.draw_networkx_edge_labels(
    G, 
    positions, 
    edge_labels={(u, v): d['label'] for u, v, d in G.edges(data=True)}
)

plt.savefig('figures/figure_1.png')