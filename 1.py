## @file 1.py
#  @brief This file contains the implementation of the RNA folding problem using dynamic programming.
import networkx as nx
import matplotlib.pyplot as plt

def visualize_sequence(sequence, pairings):
    """! 
    This function visualizes the secondary structure of an RNA sequence along with its pairings using NetworkX and Matplotlib.

    @param   sequence (str): The RNA sequence.
    @param   pairings (list of tuple): List of tuples representing pairings.
        Each tuple contains the indices of bases forming a pair.
    @return None
    """
    ## Define colors for different bases
    node_colors = {
        "A": "#FF6B7B",  # Red
        "C": "#FFD92D",  # Yellow
        "G": "#6BCB87",  # Green
        "U": "#4D96EF"   # Blue
    }

    ## Create an empty graph
    G = nx.Graph()

    ## Add nodes to the graph for each base in the sequence
    for i, c in enumerate(sequence):
        G.add_node(i+1, base=c, color=node_colors[c])

    ## Add edges to represent the sequence
    for i in range(1, len(sequence)):
        G.add_edge(i, i+1, color="#00092C", style="-", weight=3)

    ## Add edges to represent the pairings
    for e in pairings:
        G.add_edge(e[0], e[1], color="#FF0000", style="--", weight=3)

    ## Draw the graph using Kamada-Kawai layout
    nx.draw_kamada_kawai(
        G,
        labels=nx.get_node_attributes(G, "base"),  # Node labels based on base
        edge_color=nx.get_edge_attributes(G, "color").values(),  # Edge colors
        style=nx.get_edge_attributes(G, "style").values(),  # Edge styles
        node_color=nx.get_node_attributes(G, "color").values(),  # Node colors
        font_color="#151D3B",  # Font color
        font_weight="bold",  # Font weight
        linewidths=5  # Line widths
    )

    ## Set edge properties for the sequence edges
    plt.gca().collections[0].set_edgecolor("#00092C")
    plt.gca().collections[0].set_linewidth(2)

    ## Set plot title
    plt.title("RNA Secondary Structure")

    # Display the plot
    plt.show()

if __name__ == "__main__":
    ## Read pairings from file
    filename = "pairings.txt"
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = lines[0].strip()
        pairings = [(int(x.split()[0]), int(x.split()[1])) for x in lines[1:]]
    
    ## Visualize the RNA secondary structure
    visualize_sequence(sequence, pairings)
