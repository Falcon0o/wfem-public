file(GLOB WFEM_BASE_HEADER
    base/assert.h
    base/call_stack.h
)

add_subdirectory(base)
add_subdirectory(parallel)
add_subdirectory(fea)
add_subdirectory(utils)


file(GLOB WFEM_UTILS_HDF5_SRC
    utils/hdf5/attribute.cpp
    utils/hdf5/data_set.cpp
    utils/hdf5/data_space.cpp
    utils/hdf5/data_type.cpp
    utils/hdf5/file.cpp
    utils/hdf5/group.cpp
    utils/hdf5/identifier.cpp
    utils/hdf5/location.cpp
    utils/hdf5/object.cpp
    utils/hdf5/property_list.cpp
)

file(GLOB WFEM_UTILS_MNF_SRC
    utils/mnf/parser.cpp
    utils/mnf/ifstream.cpp
)


file(GLOB WFEM_CONTACT_SRC
    contact/bvh_obb.cpp
    contact/quantized_bvh.cpp
    contact/optimized_bvh.cpp
)
# add_subdirectory(multibody)
add_subdirectory(numerics)
add_subdirectory(geom)
add_subdirectory(parsers)
add_subdirectory(topo)
add_subdirectory(simulink)
add_subdirectory(graph)
add_subdirectory(mbse)

set(WFEM_SOURCE_FILES 
    ${WFEM_BASE_SRC}
    ${WFEM_PARALLEL_SRC}
    ${WFEM_MBSE_SRC}
    ${WFEM_ELEMENTS_SRC}
    ${WFEM_GEOM_SRC}
    ${WFEM_NUMERICS_SRC}
    ${WFEM_GRAPH_SRC}
    ${WFEM_UTILS_SRC}
    # ${WFEM_UTILS_HDF5_SRC}
    # ${WFEM_MARKERS_SRC}
    ${WFEM_FEA_SRC}
    ${WFEM_SIMULINK_SRC}
    # ${WFEM_CONTACT_SRC}
    ${WFEM_PARSERS_SRC}
    ${WFEM_TOPO_SRC}
    PARENT_SCOPE)


# set(WFEM_SOURCE_FILES ${WFEM_SOURCE_FILES} ${WFEM_BASE_SRC} PARENT_SCOPE)
# set(WFEM_SOURCE_FILES ${WFEM_SOURCE_FILES} ${WFEM_NUMERICS_SRC} PARENT_SCOPE)
# set(WFEM_SOURCE_FILES ${WFEM_SOURCE_FILES} ${WFEM_UTILS_SRC} PARENT_SCOPE)
# set(WFEM_SOURCE_FILES ${WFEM_SOURCE_FILES} ${WFEM_MULTIBODY_SRC} PARENT_SCOPE)
