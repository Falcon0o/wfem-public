#pragma once

#include "parallel/cuda/runtime.hpp"

#include <stdint.h>

namespace wfem {
namespace geom {

class Polyhedron;


class LinearBVH
{
public:
    LinearBVH();
    ~LinearBVH();

    void set_stream(cudaStream_t stream) {
        m_stream = stream;
    }

    void load_faces(const Polyhedron& polyhedron);
    void build_hierarchy();

    struct Aabb {
        
        // Aabb() = default;

        float4 min;
        float4 max;

        __device__ Aabb(const float4& pnt)
        {
            min.x = max.x = pnt.x;
            min.y = max.y = pnt.y;
            min.z = max.z = pnt.z;
        }

        __device__ void merge(const float4& pnt)
        {
            min.x = pnt.x < min.x ? pnt.x : min.x;
            min.y = pnt.y < min.y ? pnt.y : min.y;
            min.z = pnt.z < min.z ? pnt.z : min.z;
            max.x = pnt.x > max.x ? pnt.x : max.x;
            max.y = pnt.y > max.y ? pnt.y : max.y;
            max.z = pnt.z > max.z ? pnt.z : max.z;
        }

        __device__ void merge(const Aabb& other)
        {
            this->merge(other.min);
            this->merge(other.max);
        }
    };

    void to_obj_file(const char* file_name, int max_depth = 0) const;
private:

    void update_leaf_aabbs();
    void merge_aabbs();
    void assign_morton_codes();
    void radix_sort();

    void expand_morton_codes();
    void construct_internal_nodes();
    void update_internal_aabbs();


    int num_leaf_nodes() const {
        return m_num_leaf_nodes;
    }

    int num_internal_nodes() const {
        return m_num_leaf_nodes - 1;
    }

    Aabb* leaf_aabbs() {
        return m_aabbs + this->num_internal_nodes();
    }

private:

    int m_num_leaf_nodes;
    
    int* m_flags;
    uint32_t* m_morton_codes;
    uint64_t* m_morton_codes_64;

    uint32_t* m_payloads;
    uint32_t* m_parents;
    uint32_t* m_left_children;
    uint32_t* m_right_children;

    Aabb*   m_aabbs;
    float4* m_face_vertices;
    float4* m_face_vertices_h;

    cudaStream_t m_stream;
};


} // namespace geom
} // namespace wfem
