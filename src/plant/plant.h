#pragma once

#include <set>

#include <Eigen/Dense>

struct Vertex {
	int index;
    Eigen::Vector3d position;
	float radius;
    bool on_leaf;
	bool fixed;

	std::set<int> edges;
};

struct Edge {
	int index;
	std::array<int, 2> vertices;
};

struct Plant {
	int vertex_count;
	int edge_count;
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
};

