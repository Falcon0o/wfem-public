#pragma once

#include <memory>

namespace wfem {
namespace fea {

class Simulation
{
public:

    Simulation();
    ~Simulation();

    struct Parameters {
        bool USE_MPI = false;
    };

    Simulation& initialize(const Parameters& parameters, int* argc, char*** argv);
    Simulation& finalize();

    Simulation& load_inp_file(const char* inp_file);
    Simulation& run();
    
private:
    class Impl;
    Impl* m_impl;
};
} // namespace fea
} // namespace wfem
