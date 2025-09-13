file(GLOB WFEM_GEOM_HDR
    polyhedron.hpp
    intersect.hpp
)

file(GLOB WFEM_GEOM_SRC
    linear_bvh.cpp
    linear_bvh.cu

    polyhedron.cpp

    intersect.cpp
    optimized_bvh.cpp
)

add_subdirectory(elements)

set(WFEM_GEOM_SRC 
    ${WFEM_GEOM_SRC} ${WFEM_GEOM_ELEMENTS_SRC}
    PARENT_SCOPE)

install(FILES ${WFEM_GEOM_HDR} DESTINATION ${WFEM_HEARDER_INSTALL_PREFIX}/geom)



