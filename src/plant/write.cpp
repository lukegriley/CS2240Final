#include "write.h"

#include <fstream>

void write_plant(const std::string &plant_path, const Plant &plant) {
    std::ofstream plant_ostream(plant_path);
    if (!plant_ostream.good()) {
        throw std::runtime_error("Unable to open " + plant_path);
    }

    plant_ostream << "verts " << plant.vertex_count << std::endl;
    for (const Vertex &vertex : plant.vertices) {
        // Output plant vertices
        // Convert to Y-up
        plant_ostream << vertex.index << ","
                      << vertex.position[0] << ","
                      << vertex.position[2] << ","
                      << -vertex.position[1] << ","
                      << vertex.radius << ","
                      << static_cast<int>(vertex.on_leaf) << ","
                      << static_cast<int>(vertex.fixed) << std::endl;
    }

    plant_ostream << "edges " << plant.edge_count << std::endl;
    for (const Edge &edge : plant.edges) {
        // Output plant edges
        plant_ostream << edge.vertices[0] << ","
                      << edge.vertices[1] << std::endl;
    }

    // Possibly unnecessary
    plant_ostream.close();
}

void Writer::write_S_decomp(const std::string &S_decomp_path, const Plant &plant) {
    std::ofstream S_ostream(S_decomp_path);
    if (!S_ostream.good()) {
        throw std::runtime_error("Unable to open " + S_decomp_path);
    }

    S_ostream << "phi " << plant.phi.rows() <<" " << plant.phi.cols() << std::endl;
    for (int i = 0; i < plant.phi.rows(); ++i) {
        for (int j = 0; j < plant.phi.cols(); ++j) {
            S_ostream << plant.phi(i, j) << " ";
        }
        S_ostream << std::endl;
    }
    S_ostream << "lambda " << plant.lambda.rows() << std::endl;
    for (int i = 0; i < plant.lambda.rows(); ++i) {
        S_ostream << plant.phi(i, i);
        S_ostream << std::endl;
    }

    S_ostream.close();
}

