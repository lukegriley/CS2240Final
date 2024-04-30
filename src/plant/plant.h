#pragma once

#include <Eigen/Dense>

namespace plant {

struct Plant;

struct Vertex {
    int index = -1;
    Eigen::Vector3d tail_position;
    double radius;
    bool on_leaf = false;
    bool fixed = false;

    std::vector<int> edges;
    int parent = -1;
    std::vector<int> children;

    Eigen::Vector3d direction(const Plant &plant) const;
    Eigen::Vector3d head_position(const Plant &plant) const;
};

struct Edge {
    int index = -1;
    std::array<int, 2> vertices;
};

struct Plant {
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
};

}
