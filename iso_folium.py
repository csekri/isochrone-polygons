import folium
import shapely.wkt
from osmnx import plot

with open('out.txt') as f:
    lines = [line.rstrip() for line in f]

radii = list(map(int,lines[0].split()))
del lines[0]

gpolys = []
for i in range(len(radii)):
    wkt = shapely.wkt.loads(lines[i])
    gpolys.append(wkt)

m = folium.Map(location=[46.25970421923725, 20.167449116706845], tiles='OpenStreetMap', zoom_start=10)

iso_colours = plot.get_colors(n=len(radii)+1, cmap='hot',  start=0, return_hex=True)
del(iso_colours[0])
iso_colours = list(map(lambda x: '#'+x[1]+x[2]+x[3]+x[4]+x[5]+x[6], iso_colours))
for gpoly, colour, radius in zip(reversed(gpolys), iso_colours, reversed(radii)):
    fol = folium.GeoJson(
        gpoly,
        name='radius: '+str(radius)+'min',
        style_function=lambda x, colour=colour: {'fillColor': colour, 'fillOpacity': 0.4, 'color': colour},
        highlight_function=lambda feature: {"fillcolor": "green", "color": "green"}
    )
    fol.add_to(m)


folium.TileLayer('stamentoner').add_to(m)
folium.TileLayer('stamenterrain').add_to(m)
folium.TileLayer('stamenwatercolor').add_to(m)
folium.TileLayer('cartodbdark_matter').add_to(m)
folium.LayerControl().add_to(m)

m.save('index.html')
