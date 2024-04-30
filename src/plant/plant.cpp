#include "plant/plant.h"

using namespace Eigen;

namespace plant {

Vector3d Vertex::direction(const Plant &plant) const {
    return this->tail_position - this->head_position(plant);
}

Vector3d Vertex::head_position(const Plant &plant) const {
    if (this->parent == -1) {
        return Vector3d::Zero();
    } else {
        return plant.vertices[this->parent].tail_position;
    }
}

}
