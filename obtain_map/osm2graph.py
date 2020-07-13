import osmnx as ox

G = ox.graph_from_xml('extract.osm', bidirectional=True)
ox.save_graphml(G, 'out.graphml')
