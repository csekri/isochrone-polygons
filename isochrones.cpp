#include <fstream>
#include <string>
#include <cmath>

// BOOST GRAPH LIBRARY (BGL)
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// BOOST GEOMETRY
#include <boost/geometry.hpp>

// BOOST POLYGON
#include <boost/polygon/polygon.hpp>

// BOOST PROGRAM OPTIONS
#include<boost/program_options.hpp>

// FLAGS TO SET BEHAVIOUR
bool is_log;
double simplify_quantity;
bool fill_inside;
double travel_speed;
std::vector<double> radii;
double buffer_distance;
double additional_buffer_distance;
const double default_buffer_distance{ 0.0002 };
int points_per_circle;
bool subtract_smaller;
std::string in_file;
std::string out_file;
std::string centres;
bool predecessor;


// OTHER CONSTANTS
const double MULTIPLIER{ 100000000.0 };
// radius of Earth in km
const double R{ 6371.0 };
// factor to change degrees to radians
const double TO_RAD{ 3.1415926536 / 180.0 };


/*
ATTRIBUTES AT EACH VERTEX
We could choose to include osmid and highway, but for the sake of
saving memory we do not include them as they are not needed in the
scope of the project. We leave them in comment.
*/
struct VertexProperty
{
    // longitude at the vertex
    double x;
    // latitude at the vertex
    double y;
    // std::string osmid;
    // std::string highway;

    // predecessor vertex id of the vertex after shortest distance calculation
    int pred;
    // distance from the closest source vertex
    double dist;
};

/*
ATTRIBUTES OF EACH EDGE
Similarly we ignore osmid, name, highway, maxspeed, oneway.
*/
struct EdgeProperty
{
    // std::string osmid;
    // std::string name;
    // std::string highway;
    // std::string maxspeed;
    // std::string oneway;

    // length of the edge in meters
    double length;
    // Linestring that describes all points an edge composed of
    std::string geometry;
    // the weight in our calculations, the time it takes to travel the edge with a given speed
    double time;
    // true if a the edge is within <radius> distance from the closest cource
    std::vector<bool> included;
};


// the graph type we use with directed edges
using Graph_t = boost::adjacency_list<boost::vecS,
                                    boost::vecS,
                                    boost::directedS,
                                    VertexProperty,
                                    EdgeProperty>;

using Vertex_t = Graph_t::vertex_descriptor;
using Edge_t = Graph_t::edge_descriptor;

namespace bg = boost::geometry;
namespace gtl = boost::polygon;
namespace po = boost::program_options;

// types we define for the Boost Geometry Library
typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::polygon<point_t, false> polygon_t;
typedef bg::model::multi_polygon<polygon_t> mpolygon_t;
typedef bg::model::multi_point<point_t> mpoint_t;
typedef bg::model::linestring<point_t> linestring_t;
typedef bg::model::multi_linestring<linestring_t> mlinestring_t;
typedef bg::ring_type<polygon_t>::type ring_t;

// types we define for the Boost Polygon Library
typedef gtl::polygon_with_holes_data<long int> Polygon_Holes;
typedef gtl::polygon_data<long int> Polygon_NoHoles;
typedef gtl::point_data<long int> Point;
typedef gtl::polygon_set_data<long int> PolygonSet;
template <typename T>
struct lookup_polygon_set_type { typedef PolygonSet type; };
typedef typename lookup_polygon_set_type<gtl::property_merge<long int, long int>>::type polygon_set_type;
typedef std::map<std::set<long int>, polygon_set_type> property_merge_result_type;


// my_visitor ensures that we stop when up to radius we covered the graph,
// and do not compute the distance and predicate to all vertices
struct my_visitor : boost::default_dijkstra_visitor
{
    using base = boost::default_dijkstra_visitor;
    struct done{};

    my_visitor(double dist) : distance(dist) {}

    void finish_vertex(Vertex_t v, Graph_t const& g)
    {
        if (g[v].dist > distance)
            throw done{};

        base::finish_vertex(v, g);
    }

  private:
    double distance;
};


// spherical distance between two vertices given their azimuth angle (longitude)
// and zenith angle (latitude)
double haversine_dist(double ph1, double th1, double ph2, double th2)
{
	double dx, dy, dz;
	ph1 -= ph2;
	ph1 *= TO_RAD, th1 *= TO_RAD, th2 *= TO_RAD;

	dz = sin(th1) - sin(th2);
	dx = cos(ph1) * cos(th1) - cos(th2);
	dy = sin(ph1) * cos(th1);
	return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * R;
}


// the haversine formula computes the id of the closest node in the graph
int closest_node(double x, double y, Graph_t g)
{
    double min_dist = 1000000000.0;
    double dist = 0.0;
    int min_id;
    for (auto vd : boost::make_iterator_range(vertices(g))) {
        dist = haversine_dist(x, y, g[vd].x, g[vd].y);
        if (dist < min_dist){
            min_dist = dist;
            min_id = vd;
        }
    }
    return min_id;
}


// converts BGL type polygon into GTL type polygon
// the credit goes to
// https://github.com/mmccoo/nerd_mmccoo/blob/master/boost_polygon_geometry/poly_utils.cpp
void geom_poly2gtl_poly(const polygon_t& boost_poly, Polygon_Holes& gtlpoly)
{
    std::vector<Point> pts;
    for (point_t pt : boost_poly.outer()) {
        pts.push_back(gtl::construct<Point>((long int) (pt.get<0>()*MULTIPLIER),
                                            (long int) (pt.get<1>()*MULTIPLIER)));
    }
    gtl::set_points(gtlpoly, pts.begin(), pts.end());

    std::vector<Polygon_NoHoles> holes;
    for (ring_t r: boost_poly.inners()) {
        std::vector<Point> pts;
        for (point_t pt : r) {
            pts.push_back(gtl::construct<Point>((long int) (pt.get<0>()*MULTIPLIER),
                                                (long int) (pt.get<1>()*MULTIPLIER)));
        }
        Polygon_NoHoles hole;
        gtl::set_points(hole, pts.begin(), pts.end());
        holes.push_back(hole);
    }
    gtl::set_holes(gtlpoly, holes.begin(), holes.end());
}


// converts GTL type polygon into BGL type polygon
// the credit goes to
// https://github.com/mmccoo/nerd_mmccoo/blob/master/boost_polygon_geometry/poly_utils.cpp
void gtl_poly2geom_poly(const Polygon_Holes& gtlpoly, polygon_t& boost_poly)
{
    boost_poly.clear();
    for(Point pt : gtlpoly) {
        boost_poly.outer().push_back(
            point_t((double)gtl::x(pt)/MULTIPLIER,(double)gtl::y(pt)/MULTIPLIER));
    }

    int num_holes{ 0 };
    for(auto iter = gtlpoly.begin_holes(); iter != gtlpoly.end_holes(); ++iter) {
        num_holes++;
        const gtl::polygon_data<long int> h = *iter;

        ring_t r;
        for(Point pt : h) {
            r.push_back(point_t((double)gtl::x(pt)/MULTIPLIER,(double)gtl::y(pt)/MULTIPLIER));
        }
        boost_poly.inners().push_back(r);
    }

    // if gtl goes in the wrong direction. clockwise vs counter-clockwise.
    bg::correct(boost_poly);
}


// rescales Linestring
void scale_geom(mlinestring_t mls, mlinestring_t &mls_scaled, double scale_x, double scale_y)
{
    mls_scaled.clear();

    for (linestring_t ls : mls) {
        std::vector<point_t> pts;
        linestring_t ls_out;

        for (point_t pt : ls){
            pts.push_back(point_t(pt.get<0>()*scale_x, pt.get<1>()*scale_y));
        }
        bg::assign_points(ls_out, pts);
        mls_scaled.push_back(ls_out);
    }
}

// rescales Polygon
void scale_geom(polygon_t  poly, polygon_t &poly_scaled, double scale_x, double scale_y)
{
    poly_scaled.clear();

    for (point_t pt : poly.outer()) {
        poly_scaled.outer().push_back(point_t(pt.get<0>()*scale_x, pt.get<1>()*scale_y));
    }

    std::vector<Polygon_NoHoles> holes;
    for (ring_t r: poly.inners()) {
        std::vector<point_t> pts;
        for (point_t pt : r) {
            pts.push_back(point_t(pt.get<0>()*scale_x, pt.get<1>()*scale_y));
        }
        ring_t ring;
        bg::assign_points(ring, pts);
        poly_scaled.inners().push_back(ring);
    }
}

// rescales Multipolygon
void scale_geom(mpolygon_t  mpoly, mpolygon_t &mpoly_scaled, double scale_x, double scale_y)
{
    mpoly_scaled.clear();

    for (polygon_t poly : mpoly) {
        polygon_t poly_scaled;
        scale_geom(poly, poly_scaled, scale_x, scale_y);
        mpoly_scaled.push_back(poly_scaled);
    }
}


// reads graphml file into Graph_t
// Again we omit certain OSM properties to save memory. These are left in comment.
Graph_t ReadGraph(std::string fn)
{
    std::ifstream is(fn.c_str());
    if (!is.is_open()) {
        std::cout << "loading file '" << fn << "'failed." << std::endl;
        throw "Could not load file.";
    }
    Graph_t graph;
    boost::dynamic_properties dp(boost::ignore_other_properties);
    dp.property("x", boost::get(&VertexProperty::x, graph));
    dp.property("y", boost::get(&VertexProperty::y, graph));
    // dp.property("osmid", boost::get(&VertexProperty::osmid, graph));
    // dp.property("highway", boost::get(&VertexProperty::highway, graph));

    // dp.property("osmid", boost::get(&EdgeProperty::osmid, graph));
    // dp.property("name", boost::get(&EdgeProperty::name, graph));
    // dp.property("highwaye", boost::get(&EdgeProperty::highway, graph));
    // dp.property("maxspeed", boost::get(&EdgeProperty::maxspeed, graph));
    // dp.property("oneway", boost::get(&EdgeProperty::oneway, graph));
    dp.property("length", boost::get(&EdgeProperty::length, graph));
    dp.property("geometry", boost::get(&EdgeProperty::geometry, graph));

    boost::read_graphml(is, graph, dp);

    return graph;
};


// read the source points from text file which contains a WKT Linestring
linestring_t read_src_points(std::string centres)
{
    std::ifstream input(centres);
    std::string input_line;
    for(std::string line; getline(input, line); ) {
        input_line += line;
    }
    linestring_t src_points;
    bg::read_wkt(input_line, src_points);
    return src_points;
}


// writes radii and Multipolygons into files
void write_isochrones(std::string out_file, std::vector<mpolygon_t> isochrone_polys)
{
    std::ofstream poly_file(out_file);
    for (int r : radii) {
        poly_file << r << " ";
    }
    poly_file << std::endl;
    for (int i=0; i<radii.size(); i++) {
        poly_file << std::setprecision(12) << bg::wkt(isochrone_polys[i]) << std::endl;
    }
      poly_file.flush();
      poly_file.close();
}


// only changes the "included" property if the edge
// is reached within a certain radius with, this is done with Dijkstra's algorithm
void get_isochrones(int source, std::vector<double> radii, Graph_t &g) {
    double maxi = -1;
    for (auto radius : radii) {
        if (maxi < radius) maxi = radius;
    }
    my_visitor vis { maxi };

    try {
    boost::dijkstra_shortest_paths(g, source,
                          boost::visitor(vis)
                          .predecessor_map(get(&VertexProperty::pred, g))
                          .distance_map(get(&VertexProperty::dist, g))
                          .weight_map(get(&EdgeProperty::time, g)));
    } catch(my_visitor::done const&) {
    }


    for (int i=0; i<radii.size(); i++) {
        Vertex_t src, tar;
        auto es = boost::edges(g);
        for (auto eit = es.first; eit != es.second; ++eit) {
            src = boost::source(*eit, g);
            tar = boost::target(*eit, g);
            if (g[tar].dist < radii[i] && g[src].dist < radii[i]) {
                if(predecessor == true && g[tar].pred != src) {
                    continue;
                }
                g[*eit].included[i] = true;
            }
        }
    }
}


// unary union from a property_merge into mpolygon_t
void polygon_merge(gtl::property_merge<long int, long int> pm, mpolygon_t &multipolygon)
{
    property_merge_result_type merge_result;
    pm.merge(merge_result);
    PolygonSet result = merge_result.begin()->second;
    std::vector<Polygon_Holes> polys;
    result.get(polys);
    for (int k=0; k<polys.size(); k++) {
        polygon_t poly;
        gtl_poly2geom_poly(polys[k], poly);
        multipolygon.push_back(poly);
    }
}


// computes the Multilinestring within a given radius
void get_multilinestring_within(int i, Graph_t g, mlinestring_t &mls)
{
    auto es = boost::edges(g);

    // for a directed edge src is the source and tar is the target vertex
    Vertex_t src, tar;

    linestring_t ls;

    // we loop through all the edges in the graph
    for (auto eit = es.first; eit != es.second; ++eit) {
        // if the i-th smallest radius is in reach
        if (g[*eit].included[i] == true) {
            // if i > 0 and we already processed the edge in the previous run
            // of the outmost loop then we continue in the loop
            if (i > 0 && g[*eit].included[i-1] == true) {
                continue;
            }

            src = boost::source(*eit, g);
            tar = boost::target(*eit, g);
            std::string str(g[*eit].geometry);

            // if the geometry of the edge is not empty then we read
            // the geometry into Linestring
            if (str.length() > 1) {
                bg::read_wkt(str, ls);
            }
            // else we consider the endpoints of of the edge as Linestring
            else {
                std::vector<point_t> pts;
                pts.push_back(point_t(g[src].x, g[src].y));
                pts.push_back(point_t(g[tar].x, g[tar].y));
                bg::assign_points(ls, pts);
            }
            mls.push_back(ls);
        }
    }
}


// applies transformations on the geometry
void build_isochrone_geometry(std::vector<mpolygon_t> &isochrone_polys,
                            std::vector<Polygon_Holes> polys, double lat_lon_ratio)
{
    bg::strategy::buffer::distance_symmetric<double> distance_strategy(buffer_distance);
    bg::strategy::buffer::join_round join_strategy(points_per_circle);
    bg::strategy::buffer::end_round end_strategy(points_per_circle);
    bg::strategy::buffer::point_circle circle_strategy(points_per_circle);
    bg::strategy::buffer::side_straight side_strategy;

    // auxiliary variables for further computations with polygons
    polygon_t    poly, poly_aux;
    ring_t       ring;
    mpolygon_t   mpoly, mpoly_scaled, mpoly_collect;

    // we loop over radii starting with the smallest
    for (int j=0; j<polys.size(); j++) {
        gtl_poly2geom_poly(polys[j], poly);
            if (fill_inside == true && predecessor == false) {
                ring = poly.outer();
                bg::convert(ring, poly_aux);
                bg::simplify(poly_aux, poly, simplify_quantity);
                bg::strategy::buffer::distance_symmetric<double>
                    distance_strategy(additional_buffer_distance);
                bg::buffer(poly, mpoly, distance_strategy, side_strategy,
                           join_strategy, end_strategy, circle_strategy);
                scale_geom(mpoly, mpoly_scaled, 1/lat_lon_ratio, 1);
                ring.clear();
                ring = mpoly_scaled[0].outer();
                mpoly_scaled.clear();
                mpoly.clear();
                poly.clear();
                bg::convert(ring, poly);
            } else {
                bg::simplify(poly, poly_aux, simplify_quantity);
                scale_geom(poly_aux, poly, 1/lat_lon_ratio, 1);
            }
            mpoly_collect.push_back(poly);
    }
    isochrone_polys.push_back(mpoly_collect);
}


void union_of_fill_inside_true_predecessor_false(std::vector<mpolygon_t> &isochrone_polys)
{
    std::vector<mpolygon_t> isochrone_polys_union;
    gtl::property_merge<long int, long int> pm;
    for (mpolygon_t mpoly : isochrone_polys) {
        for (polygon_t poly : mpoly){
            Polygon_Holes pl;
            geom_poly2gtl_poly(poly, pl);
            pm.insert(pl, 0);
        }
        mpolygon_t multipolygon, multipolygon_outer;
        polygon_merge(pm, multipolygon);
        for (polygon_t poly: multipolygon) {
            ring_t ring = poly.outer();
            poly.clear();
            bg::convert(ring, poly);
            multipolygon_outer.push_back(poly);
        }
        isochrone_polys_union.push_back(multipolygon_outer);
    }
    isochrone_polys.clear();
    isochrone_polys.insert(isochrone_polys.end(),
                           isochrone_polys_union.begin(), isochrone_polys_union.end());
}


void subtract_smaller_polygon(std::vector<mpolygon_t> &isochrone_polys)
{
    std::vector<mpolygon_t> isochrone_polys_difference;
    for (int j=0; j<isochrone_polys.size()-1; j++) {
        if (j == 0) {
            isochrone_polys_difference.push_back(isochrone_polys[0]);
        }
        gtl::property_merge<long int, long int> pm;
        Polygon_Holes pl;
        for (polygon_t poly : isochrone_polys[j+1]) {
            geom_poly2gtl_poly(poly, pl);
            pm.insert(pl, 0);
        }
        for (polygon_t poly : isochrone_polys[j]) {
            geom_poly2gtl_poly(poly, pl);
            pm.insert(pl, 1);
        }
        mpolygon_t multipolygon;
        polygon_merge(pm, multipolygon);
        isochrone_polys_difference.push_back(multipolygon);
    }
    isochrone_polys.clear();
    isochrone_polys.insert(isochrone_polys.end(),
                           isochrone_polys_difference.begin(), isochrone_polys_difference.end());
}


// obtaining command line options
int cmd_input(int ac, char * av[],
            bool &is_log,
            double &simplify_quantity,
            bool &fill_inside,
            double &travel_speed,
            std::vector<double> &radii,
            double &buffer_distance,
            double &additional_buffer_distance,
            int &points_per_circle,
            bool &subtract_smaller,
            std::string &in_file,
            std::string &out_file,
            std::string &centres,
            bool &predecessor)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("verbose", po::bool_switch(&is_log)->default_value(false),
        "prints log details")
        ("simplify", po::value<double>(&simplify_quantity),
        "the epsilon in Ramer-Douglas-Peucker algorithm, the bigger the number is, the less points shall be present in the output")
        ("fill", po::bool_switch(&fill_inside)->default_value(false),
        "if the polygon contains holes it removes them")
        ("travel-speed", po::value<double>(&travel_speed)->default_value(20.0),
        "the travel speed in km/h")
        ("radii", po::value<std::vector<double>>(&radii)
                ->multitoken()->default_value(std::vector<double>{10,20,30,40,50,60},
        "10,20,30,40,50,60"), "the radii of the isochrone curves in minutes")
        ("circle-points", po::value<int>(&points_per_circle)->default_value(1),
        "the radii of the isochrone curves in minutes")
        ("subtract", po::bool_switch(&subtract_smaller)->default_value(false),
        "subtracts the smaller radius polygon from the larger thus creating a hole in that")
        ("input-file", po::value<std::string>(&in_file)->required(),
        "the input .graphml file")
        ("output-file,o", po::value<std::string>(&out_file)->default_value("out.txt"),
        "the output text file")
        ("sources", po::value<std::string>(&centres)->default_value("sources.txt"),
        "a file containing a WKT Linestring whose vertices are the sources")
        ("predecessor", po::bool_switch(&predecessor)->default_value(false),
        "computes the isochrone curves with a shortest path tree, fill does not affect it, it is not adviced to have multiple sources in this case")
        ("buffer-distance", po::value<double>(&buffer_distance)->default_value(default_buffer_distance),
        "thickness added to the roads and polygons")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    po::notify(vm);

    if (buffer_distance <= 0) {
        std::cout << "the buffer distance must be positive\n";
        return 1;
    }

    if (fill_inside == true and predecessor == false) {
        if (buffer_distance > default_buffer_distance*2) {
            additional_buffer_distance = buffer_distance - default_buffer_distance;
            buffer_distance = default_buffer_distance;
        }
        else {
            buffer_distance /= 2;
            additional_buffer_distance = buffer_distance;
        }
    }
    else {
        if (buffer_distance > default_buffer_distance*2) {
            std::cout << "WARNING: too large buffer distance may result in the program never complete.\n";
        }
    }
    return 0;
}


int main(int ac, char * av[])
{
    // command line arguments
    int error = cmd_input(
            ac, av,
            is_log,
            simplify_quantity,
            fill_inside,
            travel_speed,
            radii,
            buffer_distance,
            additional_buffer_distance,
            points_per_circle,
            subtract_smaller,
            in_file,
            out_file,
            centres,
            predecessor
    );

    if (error == 1){
        return 1;
    }

    if (is_log == true) {
        std::cout << "File \"" << in_file << "\" starts being read.\n";
    }

    Graph_t g = ReadGraph(in_file);

    if (is_log == true) {
        std::cout << "File \"" << in_file << "\" is read.\n";
    }

    double meters_per_minute{ travel_speed * 1000.0 / 60.0 };
    auto es = boost::edges(g);

    // we initialise "included" as a vector of falses
    for (auto eit = es.first; eit != es.second; ++eit) {
        g[*eit].time =  (g[*eit].length / meters_per_minute);
        for (int i=0; i<radii.size(); i++) {
            g[*eit].included.push_back(false);
        }
    }

    // 40075==circumference of Earth in km;
    // 111.32==1 degree of latitude in km
    // lat_lon_ratio expresses the ratio of one meter of latitude and longitude
    // at given latitude
    double lat_lon_ratio{ 40075.0 * cos(g[0].y * TO_RAD) / 360.0 / 111.32 };

    // loads the source points
    linestring_t src_points = read_src_points(centres);

    int idx = 0;
    auto size = src_points.size();
    for (point_t pt : src_points) {
        idx++;
        int source = closest_node(pt.get<0>(), pt.get<1>(), g);
        if (is_log == true) {
            std::cout << "#" <<  idx << "/" << size << ", " << "x="
                      << pt.get<0>() << ", y=" << pt.get<1>() << "\n";
        }
        get_isochrones(source, radii, g);
    }
    if (is_log == true) {
        std::cout << "Dijkstra is done.\n";
        std::cout << "Geometry is being assembled.\n";
    }

    //buffer strategy, for reference visit the website of BGL
    bg::strategy::buffer::distance_symmetric<double> distance_strategy(buffer_distance);
    bg::strategy::buffer::join_round join_strategy(points_per_circle);
    bg::strategy::buffer::end_round end_strategy(points_per_circle);
    bg::strategy::buffer::point_circle circle_strategy(points_per_circle);
    bg::strategy::buffer::side_straight side_strategy;


    // we collect the computed isochronic polygons in "isochrone_polys"
    std::vector<mpolygon_t> isochrone_polys;

    // the result of GTL property merge (union in our case) is converted into
    // std::vector<Polygon_Holes> which is stored in "polys,"
    // futhermore the previous (smaller) isochronic curve is inserted into the
    // next one thus speeding up computation, which is done with "polys"
    std::vector<Polygon_Holes> polys;

    // we loop over radii starting with the smallest
    for (int i=0; i<radii.size(); i++) {
        // auxiliary variables for polygon manipulations
        mlinestring_t mls, mls_scaled;
        mpolygon_t   mpoly_aux;
        Polygon_Holes pl;

        get_multilinestring_within(i, g, mls);

        // at this point if collected in a vector, the isochronic Multilinestrings
        // could be saved

        // we collect the polygons to be unioned in "pm"
        gtl::property_merge<long int, long int> pm;

        // we insert the previous isochronic curve here thus gaining speed
        for (Polygon_Holes polygon : polys) {
            pm.insert(polygon, 0);
        }
        // polys emptied for further computations after insertion
        polys.clear();

        // in buffer we scale the geometries thus buffer affects x and y directions equally
        scale_geom(mls, mls_scaled, lat_lon_ratio, 1);

        // memory saving for the rest of the loop
        mls.clear();

        // buffer applied to the scaled Multilinestring
        for (linestring_t ls : mls_scaled) {
            bg::buffer(ls, mpoly_aux, distance_strategy, side_strategy,
                       join_strategy, end_strategy, circle_strategy);
            geom_poly2gtl_poly(mpoly_aux.at(0), pl);
            pm.insert(pl, 0);
        }

        // "merge_result" holds the result for the unary union
        property_merge_result_type merge_result;
        pm.merge(merge_result);
        PolygonSet result = merge_result.begin()->second;
        result.get(polys);

        build_isochrone_geometry(isochrone_polys, polys, lat_lon_ratio);

    }

    if (fill_inside == true && predecessor == false) {
        union_of_fill_inside_true_predecessor_false(isochrone_polys);
    }

    if (subtract_smaller == true) {
        subtract_smaller_polygon(isochrone_polys);
    }

    // writing result into file
    write_isochrones(out_file, isochrone_polys);

      if (is_log == true) {
          std::cout << "Polygons are written into \"" << out_file << "\"\n";
      }
      return 0;
}
