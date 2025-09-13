#include <iostream>
// #include <dlfcn.h>



#include "fea/simulation.hpp"

// #include <mpi.h>

int main(int argc, char** argv)
{
    wfem::fea::Simulation sim;

    wfem::fea::Simulation::Parameters params;
    params.USE_MPI = false;

    sim.initialize(params, &argc, &argv);
    sim.load_inp_file("/home/wangxinyu/workspace/myself/wfem/data/conrod.inp");
    sim.run();

    sim.finalize();

    return 0;
}

    // constexpr bool USE_MPI = false;

    // if (USE_MPI) { wfem::MPI::init(&argc, &argv); }

    // wfem::Logger::initialize("wfem");

    // constexpr bool single_precision = true;
    // auto physical_model = wfem::fea::PhysicalModel::create(single_precision);

    // physical_model->

    // if (USE_MPI) { wfem::MPI::finalize(); };