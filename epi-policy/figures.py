import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

############################################################
# Figure 1
############################################################

# Use default matplotlib mathtext renderer instead of LaTeX
plt.rcParams['text.usetex'] = False


nodes = ["S", "E", "P", "I", "A", "R", "D"]

# Instantiate nodes for directed graph
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

############################################################
# Figure 2
############################################################

results_long = pd.read_csv("/Users/abhay/Documents/XLab/epi-policy/results/raw_results_long.csv")
MA_deaths = pd.read_csv("/Users/abhay/Documents/XLab/epi-policy/data/data_table_for_cumulative_deaths__massachusetts.csv", header = 2)

MA_deaths['Date'] = pd.to_datetime(MA_deaths['Date'], format='%b %d %Y')
MA_deaths['Cumulative Deaths'] = MA_deaths['Cumulative Deaths'].replace('Counts 1-9', '0').astype(int)
MA_deaths['Days'] = (MA_deaths['Date'] - MA_deaths['Date'].min()).dt.days
MA_deaths = MA_deaths[MA_deaths['Days From Start'] <= 365]

# Existing code to generate a plot (assuming df is defined elsewhere)
total_deaths = results_long.groupby('day')['D'].sum()

plt.figure(figsize=(10, 6))
plt.plot(total_deaths.index, total_deaths.values, marker='o', linestyle='-', label='Total Deaths by Day')

# Plot the 'Cumulative Deaths' from df_new
plt.plot(MA_deaths['Days'], MA_deaths['Cumulative Deaths'], marker='x', linestyle='--', label='Deaths')

# Enhancing the plot
plt.title('Comparison of Total Deaths Over Time')
plt.xlabel('Days')
plt.ylabel('Number of Deaths')
plt.legend()
plt.grid(True)
plt.show()