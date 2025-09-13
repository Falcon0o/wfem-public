#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <string_view>
#include "Eigen/Dense"

namespace wfem {
namespace geom {


/// @brief 三维表面网格，包含拓扑关系
class Polyhedron 
{
public:
    Polyhedron(const Polyhedron&);

    class HalfEdge;

    class Vertex 
    {
    public:

        Vertex() :
            m_point(nullptr),
            m_half_edge(nullptr)
        {}

        HalfEdge* half_edge() const { return m_half_edge; }
        void set_half_edge(HalfEdge* he) { m_half_edge = he; }
        int degree() const;

        const Eigen::Vector3d& point() const { return *reinterpret_cast<const Eigen::Vector3d*>(m_point); }
        void set_point(const double* pnt) { m_point = pnt; }

        void for_each_half_edge(std::function<void(const HalfEdge*)> callback) const;

    private:

        HalfEdge* m_half_edge;
        const double* m_point;
    }; // class Vertex

    class Face 
    {
    public:
        Face() :
            m_half_edge(nullptr)
        {}

        Face(const Face& other) 
        :   m_half_edge(other.m_half_edge) 
        { }

        const HalfEdge* half_edge() const { return m_half_edge; }
        void set_half_edge(const HalfEdge* he) { m_half_edge = he; }
        
        int degree() const;

        void for_each_half_edge(std::function<void(const HalfEdge*)> callback) const;
        // void for_each_vertex(std::function<void(Vertex*)> callback) const;

        bool is_neighbor(const Face* face) const;

        const HalfEdge* common_half_edge(const Face*) const;

    private:

        const HalfEdge* m_half_edge;
    }; // class Face

    /// @brief 
    class HalfEdge 
    {
    public:
        HalfEdge() :
            m_opposite(nullptr),
            m_prev(nullptr),
            m_next(nullptr),
            m_vertex(nullptr),
            m_face(nullptr) 
        {

        }
#if 0
        HalfEdge& operator=(const HalfEdge& other) 
        {
            if (this != &other)
            {
                m_opposite = other.m_opposite;
                m_prev = other.m_prev;
                m_next = other.m_next;
                m_vertex = other.m_vertex;
                m_face = other.m_face;
            }
            return *this;
        }
#endif
        const HalfEdge* prev() const { return m_prev; }
        void set_prev(const HalfEdge* he) { m_prev = he; }

        const HalfEdge* next() const { return m_next; }
        void set_next(const HalfEdge* he) { m_next = he; }

        const HalfEdge* opposite() const { return m_opposite; }
        void set_opposite(const HalfEdge* he) { m_opposite = he; }

        const Vertex* vertex() const { return m_vertex; }
        void set_vertex(const Vertex* vtx) { m_vertex = vtx; }

        const Face* face() const { return m_face; }
        void set_face(const Face* f) { m_face = f; }

    private:
        
        const HalfEdge* m_opposite; 
        const HalfEdge* m_prev;
        const HalfEdge* m_next;

        /// @brief 指向的点
        const Vertex* m_vertex;
        const Face* m_face;
    }; // class HalfEdge

    void write_obj_file(const char* path, bool show_face_connection = false) const;
    
    Polyhedron(const char* name);
    ~Polyhedron();

    std::size_t size_of_vertices() const { return m_size_of_vertices; }
    std::size_t size_of_half_edges() const { return m_size_of_half_edges; }
    std::size_t size_of_faces() const { return m_size_of_faces; }

    void reserve_vertices(std::size_t);
    void reserve_half_edges(std::size_t);
    void reserve_faces(std::size_t);

    Vertex* emplace_back_vertex();
    HalfEdge* emplace_back_half_edge();
    Face* emplace_back_face();

    const HalfEdge* half_edge(int i) const { return m_half_edges + i; }
    HalfEdge* half_edge(int i) { return m_half_edges + i; }

    const Vertex* vertex(int i) const { return m_vertices + i; }
    Vertex* vertex(int i) { return m_vertices + i; }

    const Face* face(int i) const { return m_faces + i; }
    Face* face(int i) { return m_faces + i; }

    int index(const Vertex* vertex) const { return vertex - m_vertices; }
    int index(const HalfEdge* half_edge) const { return half_edge - m_half_edges; }
    int index(const Face* face) const { return face - m_faces; }

    // void for_each_face(std::function<void(const Face*)> callback);

    bool has(const Face* face) const {
        return face >= m_faces & face < (m_faces + m_size_of_faces); 
    }

    bool has(const HalfEdge* edge) const {
        return edge >= m_half_edges & edge < (m_half_edges + m_size_of_half_edges); 
    }
    const std::string_view name() 
    {
        return m_name;
    }
    /// @brief 返回的半边和输入的半边都指向相同的新顶线
    HalfEdge* split_edge(HalfEdge*);
private:

private:
    std::string             m_name;

    HalfEdge*               m_half_edges;
    Vertex*                 m_vertices;
    Face*                   m_faces;

    std::size_t             m_capacity_of_half_edges;             
    std::size_t             m_capacity_of_vertices;             
    std::size_t             m_capacity_of_faces;             

    std::size_t             m_size_of_half_edges;             
    std::size_t             m_size_of_vertices;             
    std::size_t             m_size_of_faces;
};

std::vector<std::shared_ptr<Polyhedron>> read_obj_file(const char* path, std::vector<double>& vertices);
} // namespace geom
} // namespace wfem
