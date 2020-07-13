import osmnx as ox
ox.config(log_console=True)

source_coord = (46.252382, 20.148722)
G = ox.graph_from_point(source_coord, dist=100, network_type='bike')
ox.save_graphml(G, 'out.graphml')
