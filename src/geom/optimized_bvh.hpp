#pragma once


#include <vector>

#include "geom/aabb.cuh"
#include <functional>

namespace wfem {
namespace geom {

class Polyhedron;
template <typename FLOAT> struct PrimitiveBox;


//! @note 实际上，目前这里的16字节对齐是没什么卵用的，因为未使用SSE指令集

template <typename FLOAT>
class OptimizedBVH
{
public:

    OptimizedBVH();
    ~OptimizedBVH();

public:

    /// @brief OptimizedBvhNode contains both internal and leaf node information.
    ///        Total node size is 44 bytes / node. You can use the compressed version of 16 bytes.
    struct Node
    {
        Aabb<FLOAT> m_aabb;

        int m_escape_index_or_data_index;

	public:

		Node() : m_escape_index_or_data_index(0) {}

		bool is_leaf_node() const { return m_escape_index_or_data_index >= 0; }
		bool is_internal_node() const { return !this->is_leaf_node(); }

		int escape_index() const { return -m_escape_index_or_data_index; }
		void set_escape_index(int index) { m_escape_index_or_data_index = -index; }

		int data_index() const { return m_escape_index_or_data_index; }
		void set_data_index(int index) { m_escape_index_or_data_index = index; }

        bool overlap(const Node& other) const 
        {
            return m_aabb.overlap(other.m_aabb);
        }

    };
    // static_assert(sizeof(OptimizedBVH::Node) == 32);

    void build_hierarchy(const Polyhedron&);

    void to_obj_file(const char* file_name, int max_depth = 0) const;
	int depth() const { return depth(0);  }
    
    const Node& node(int index) const { return m_nodes[index]; }

    int left_node_index(int node_index) const 
    { 
        return node_index + 1; 
    }

    int right_node_index(int node_index) const
    {
        const auto& next_node = m_nodes[node_index + 1];
        if (next_node.is_leaf_node())
        {
            return node_index + 2;
        }
        return node_index + 1 + next_node.escape_index();
    }

private:
    void build_sub_tree(std::vector<PrimitiveBox<FLOAT>>& primitive_boxes, int first_index, int last_index);

    //! @brief Find Best Splitting Axis and where to split it. Sort the incoming 'leafNodes' array within range 'startIndex/endIndex'.
    int find_split_index_and_sort(std::vector<PrimitiveBox<FLOAT>>& primitive_boxes, int first_index, int last_index);

    int depth(int node_index) const;

private:

    int m_num_nodes = 0; 
    std::vector<OptimizedBVH::Node> m_nodes;
};

} // namespace geom
} // namespace wfem

template <typename FLOAT>
void find_collision(const wfem::geom::OptimizedBVH<FLOAT>& first_bvh, const wfem::geom::OptimizedBVH<FLOAT>& second_bvh, std::function<void(int, int)> handler);