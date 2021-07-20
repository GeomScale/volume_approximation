#include <iostream>
#include <fstream>
#include "misc.h"
#include "poset.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "cartesian_geom/point.h"
#include "orderpolytope.h"

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

typedef typename Poset::RT RT;
typedef typename Poset::RV RV;

Poset read_poset_from_file(std::string filename) {
    std::ifstream data_file;
    data_file.open(filename);

    unsigned int n;
    data_file >> n;

    RT curr_relation;
    RV relations;
    while(data_file >> curr_relation.first >> curr_relation.second)
        relations.push_back(curr_relation);

    return Poset(n, relations);
}


int main(int argc, char const *argv[]) {
    std::cout << "\nPoset operations: \n";
    //  ----------- basic poset operations -----------
    Poset poset = read_poset_from_file(std::string(argv[1]));
    poset.print();

    std::cout << "Checking if a sample Point (linearly spaced coordinates) lies inside the poset: ";
    std::vector<double> temp = linspace((double)(0.0), (double)(1.0), poset.num_elem());
    std::cout << poset.is_in(temp) << std::endl;
    std::cout << "\n";
    //  ----------------------------------------------


    // -------------- order polytope operations --------
    std::cout << "\nOrder polytope operations: \n";
    typedef Cartesian<double> Kernel;
    typedef typename Kernel::Point Point;
    OrderPolytope<Point> OP(poset);
    OP.print();
    std::cout << "\n";


    std::cout << "intersection of the order polytope with ray from (0.5, 0.5 .... 0.5) towards the origin" << std::endl;
    Point origin(OP.dimension());
    Point start_point(OP.dimension(), std::vector<double>(OP.dimension(), 0.5));
    Point direction = origin - start_point;
    std::pair<double, double> curr_res = OP.line_intersect(start_point, direction, true);
    Point intersect_point = start_point + curr_res.first * direction;
    intersect_point.print();
    std::cout << "\n";


    std::cout << "distances of all hyperplanes from origin: " << std::endl;
    for (const auto val: OP.get_dists(0.0))
        std::cout << val << ' ';
    std::cout << "\n\n";


    OP.normalize();
    std::cout << "normalized order polytope: " << std::endl;
    OP.print();
    std::cout << "\n";


    std::cout << "distances of all hyperplanes from origin (after normalization): " << std::endl;
    for (const auto val: OP.get_dists(0.0))
        std::cout << val << ' ';
    std::cout << "\n\n";


    std::cout << "compute reflection (requires normalization) of an incident ray with the facet number 2d (the first relation facet)" << std::endl;
    Point ray = Point::all_ones(OP.dimension());
    ray.set_coord(0, 1.5);
    std::cout << "incident ray: ";
    ray.print();
    
    OP.compute_reflection(ray, Point(), 2*OP.dimension());
    std::cout << "reflected ray: ";
    ray.print();
    // ---------------------------------------------    

    return 0;
}
