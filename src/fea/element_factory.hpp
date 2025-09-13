#pragma once

#include <vector>
#include <map>
#include <memory>

namespace wfem {
enum class AbaqusElementType;

namespace fea {

class Element;
class Elements;

class ElementFactory {
public:
    void create(std::vector<Element*>& elements, AbaqusElementType type, int n_elements, const int node_indices[]);

    template <typename ELEMENT>
    void create(std::vector<Element*>& elements, int n_elements, const int node_indices[]);

private:
    std::map<AbaqusElementType, std::shared_ptr<Elements>> m_elements;
};
} // namespace fea
} // namespace wfem
