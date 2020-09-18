# https://stackoverflow.com/questions/60147431/how-to-put-a-label-on-a-country-with-python-cartopy
# https://stackoverflow.com/questions/59276335/how-to-draw-edge-weights-using-a-weighted-adjacency-matrix
import matplotlib.pyplot as plt
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import numpy as np
import networkx as nx

def get_rgb(value):
    r = min(1, 2 * 1 * (value))
    g = min(1, 2 * 1 * (1 - value))
    b = 0
    return (r, g, b)

def plot_map(country_data, compartment, title):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.95)
    ax.add_feature(cartopy.feature.RIVERS)

    ax.set_extent([-150, 60, -25, 60])

    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    for country in countries:
        country_name = country.attributes['SOVEREIGNT']

        if country_name in country_data.keys():
            rgb_color = get_rgb(country_data[country_name])
            g = ax.add_geometries(country.geometry, ccrs.PlateCarree(), facecolor=rgb_color, label=country_name)
            x = country.geometry.centroid.x
            y = country.geometry.centroid.y
            ax.text(x, y, compartment + " = " + str(country_data[country_name]), color='black', size=15, ha='center', va='center', transform=ccrs.PlateCarree())
        else:
            ax.add_geometries(country.geometry, ccrs.PlateCarree(), facecolor=(1, 1, 1), label = country_name)

    plt.rcParams["figure.figsize"] = (50,50)
    plt.title(title)
    plt.show()

def plot_graph(W, country_data, title, scale_node_size):
    # Create DiGraph from W
    G = nx.from_numpy_matrix(W, create_using=nx.DiGraph)

    # Use spring_layout to handle positioning of graph
    layout = nx.spring_layout(G)

    # Use a list for node_sizes
    sizes = [sum(row) * scale_node_size for row in W]

    # Use a list for node colours
    color_map = []
    tot_travel = sum(sizes)
    for country_travel in sizes:
        color_map.append(get_rgb(country_travel/tot_travel))

    country_labels = {}
    country_labels.update(zip(list(range(len(country_data))), country_data.keys()))

    plt.title(title, fontsize=20)

    # Draw the graph using the layout - with_labels=True if you want node labels.
    nx.draw(G, layout,  node_size=sizes, node_color=color_map)
    nx.draw_networkx_labels(G, layout, country_labels, font_size=15)

    # Get weights of each edge and assign to labels
    labels = nx.get_edge_attributes(G, "weight")

    # Draw edge labels using layout and list of labels
    nx.draw_networkx_edge_labels(G, pos=layout, edge_labels=labels)

    # Show plot
    plt.show()

title = 'Regional levels of I'
country_data = {"Sweden": 0.8, "Norway": 0.5, "Denmark": 0.3, "Finland": 0.1}
#plot_map(country_data, "I", title)

W = np.array([[0, 2, 5, 10],
              [2, 0, 3, 3],
              [5, 3, 0, 1],
              [10, 3, 1, 0]])

title = 'Undirected Travel Graph'

#plot_graph(W, country_data, title, 1000)

