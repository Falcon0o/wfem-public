// #include "fea/elements/element_set.hpp"

#include <fstream>
#include <iostream>
#include "wfem/geom/intersect.hpp"
#include "wfem/geom/polyhedron.hpp"
// #include "base/assert.h"
#include <chrono>




int main(int argc, char** argv)
{
    std::vector<double> vertices;

    const char* obj_path = "/home/wangxinyu/workspace/myself/wfem/data/untitled.obj";
    std::vector<std::shared_ptr<wfem::geom::Polyhedron>> polyhedrons = wfem::geom::read_obj_file(obj_path, vertices);

    wfem::geom::Intersect<wfem::geom::Polyhedron, wfem::geom::Polyhedron> intersect(*polyhedrons[0], *polyhedrons[1]);
    intersect.process();

    return 0;
}