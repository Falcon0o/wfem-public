// #include "fea/elements/element_set.hpp"

#include <fstream>
#include <iostream>
#include "geom/linear_bvh.cuh"
#include "geom/polyhedron.hpp"
// #include "base/assert.h"
#include <chrono>


int main(int argc, char** argv)
{
    try {
        cudaStream_t stream = nullptr;
        WFEM_CUDA_CALL(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        {
            std::vector<double> vertices;

            const char* obj_path = "/home/wangxinyu/workspace/myself/wfem/data/untitled.obj";
            std::vector<std::shared_ptr<wfem::geom::Polyhedron>> polyhedrons = wfem::geom::read_obj_file(obj_path, vertices);

            wfem::geom::LinearBVH bvh;
            bvh.set_stream(stream);

            bvh.load_faces(*polyhedrons[1]);
            bvh.build_hierarchy();
            bvh.to_obj_file("tmp/bvh_on_gpu");
        }
        
        cudaStreamSynchronize(stream);
        cudaStreamDestroy(stream);
    }
    catch (char const* err) {
        std::cout << err << std::endl;
    }
    return 0;
}

