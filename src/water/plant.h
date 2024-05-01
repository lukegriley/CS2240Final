#pragma once

#include <set>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "plant/plant.h"

namespace water {

static float DYNAMIC_VISCOSITY = 9e-7;

struct Segment {
    int index;

    Eigen::Vector3d head_position;
    Eigen::Vector3d tail_position;
    double radius;
    bool on_leaf;
    std::vector<int> edges;

    double loss_rate;
    double volume;
    double water_amt;//theta
    double resistance;

    double height() const;
};

struct Edge {
    std::array<int, 2> segments;
    double resistance;
};

struct Plant {
    std::vector<Segment> segments;
    std::vector<Edge> edges;

    Eigen::MatrixXf phi;
    Eigen::DiagonalMatrix<float, Eigen::Dynamic> lambda;

    float delta = 3e-4; //uniform loss rate of plant

    Plant() = default;
    void initStructure(const plant::Plant &basic_plant);
    void initDiffusion();
    void updateDiffusionDelta(const float dt); // delta time diffusion update

    Eigen::MatrixXd m_omega;
    double m_alpha;
    Eigen::VectorXd m_beta;
    Eigen::SparseMatrix<double> m_S_sparse;
};


}
