#include "loader.h"

#include <fstream>

using namespace Eigen;



Plant Loader::load_plant(const std::string &plant_path) {
	std::ifstream plant_istream(plant_path);
    if (!plant_istream.good()) {
        throw std::runtime_error("Unable to open " + plant_path);
    }

	Plant plant;

	// Read the number of vertices.
    plant.vertex_count = Loader::read_header<int>(plant_istream, "verts");
	// Allocate space to hold vertex data.
	plant.vertices = std::vector<Vertex>(plant.vertex_count);

	// Read all vertex data.
	for (int i = 0; i < plant.vertex_count; ++i) {
		Vertex vertex;
        plant_istream >> vertex.index >> expect(',')
                >> vertex.tail_position[0] >> expect(',')
                >> vertex.tail_position[1] >> expect(',')
                >> vertex.tail_position[2] >> expect(',')
                >> vertex.radius >> expect(',')
                >> vertex.on_leaf;
        // Convert to z-up
        vertex.tail_position = Matrix3d {
            { 1, 0, 0 },
            { 0, 0, 1 },
            { 0, -1, 0 },
        } * vertex.tail_position;
        // Read optional arguments
        if (plant_istream.get() == ',') {
            plant_istream >> vertex.fixed;
        }

		assert(vertex.index == i);
		plant.vertices[i] = vertex;
	}

	// Read the number of edges
	plant.edge_count = read_header<int>(plant_istream, "edges");
	// Allocate space for all edges
	plant.edges = std::vector<Edge>(plant.edge_count);
	// Read all edge data.
	for (int i = 0; i < plant.edge_count; ++i) {
		Edge edge { .index = i };
        plant_istream >> edge.vertices[0] >> expect(',')
                >> edge.vertices[1];

		// Add the edge to the plant
		plant.edges[i] = edge;
		// Add the edge to all vertices
		for (int j = 0; j < 2; ++j) {
			plant.vertices[edge.vertices[j]].edges.insert(i);
		}
        // Record the first vertex's tail as the next vertex's head
        plant.vertices[edge.vertices[1]].head_position = plant.vertices[edge.vertices[0]].tail_position;
	}
    plant.initDiffusion(false);

	return plant;
}

void Loader::load_S_decomp(const std::string &plant_path, Plant &plant){
    std::ifstream plant_istream(plant_path);
    if (!plant_istream.good()) {
        throw std::runtime_error("Unable to open " + plant_path);
    }
    int phi_size = Loader::read_header<int>(plant_istream, "phi");
    plant.phi = Eigen::MatrixXf::Zero(phi_size,phi_size);

    for (int i = 0; i < phi_size; ++i) {
        for(int j=0;j<phi_size;j++) {
            float coeff;
            plant_istream >> coeff >> expect(' ');
            plant.phi(i,j) = coeff;
        }
        expect('\n');
    }

    // int lambda_size = Loader::read_header<int>(plant_istream, "lambda");
    // plant.lambda = Eigen::MatrixXf::Zero(lambda_size,lambda_size);

    // for (int i = 0; i < lambda_size; ++i) {
    //     float coeff;
    //     plant_istream >> coeff >> expect('\n');
    //     plant.lambda(i,i) = coeff;
    // }

}


