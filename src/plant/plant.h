#pragma once

#include <set>

#include <Eigen/Dense>
#include <Eigen/Sparse>

static float DYNAMIC_VISCOSITY = 1.f;

struct Vertex {
	int index;
    Eigen::Vector3d position;
    float radius;
    bool on_leaf;
	bool fixed;
	std::set<int> edges;
    float loss_rate;
    float volume;
    float water_amt;//theta
    float resistance;
};

struct Edge {
	int index;
	std::array<int, 2> vertices;
    float resistance;
};

struct Plant {
	int vertex_count;
	int edge_count;
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;

    Eigen::MatrixXf D_V;
    Eigen::MatrixXf D_l;
    Eigen::MatrixXf S;
    Eigen::MatrixXf phi;
    Eigen::MatrixXf lambda;

    float delta = 3e-4; //uniform loss rate of plant
    void initDiffusion();
    void updateDiffusion(float time);//time is total running time, not delta

    Eigen::MatrixXf theta_0;//initial water content at each node
    Eigen::MatrixXf zeta_0;//implicit computation of Phi inverse at initialization; see section 4.2
    bool load_S = true;

};





