#pragma once

#include <set>

#include <Eigen/Dense>
#include <Eigen/Sparse>

static float DYNAMIC_VISCOSITY = 1.f;

struct Vertex {
    Vertex();
    // Return the height of this cylinder
    inline double height() const {
        return (this->tail_position - this->head_position).norm();
    }

    int index;
    Eigen::Vector3d head_position;
    Eigen::Vector3d tail_position;
    double radius;
    bool on_leaf;
	bool fixed;
	std::set<int> edges;
    double loss_rate;
    double volume;
    double water_amt;//theta
    double resistance;
};

struct Edge {
	int index;
	std::array<int, 2> vertices;
    double resistance;
};

struct Plant {
	int vertex_count;
	int edge_count;
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;

    Eigen::MatrixXd D_V;
    Eigen::MatrixXd D_l;
    Eigen::MatrixXd S;
    Eigen::MatrixXf phi;
    Eigen::DiagonalMatrix<float, Eigen::Dynamic> lambda;

    float delta = 3e-4; //uniform loss rate of plant
    void initDiffusion(bool precompute = true);
    void updateDiffusion(float time);//time is total running time, not delta
    void updateDiffusionDelta(const float dt); // delta time diffusion update

    Eigen::MatrixXd m_omega;
    double m_alpha;
    Eigen::VectorXd m_beta;
    Eigen::SparseMatrix<double> m_S_sparse;

    Eigen::MatrixXf theta_0;//initial water content at each node

};



