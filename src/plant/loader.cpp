#include "loader.h"

#include <fstream>

using namespace Eigen;

// Read a line formatted as <header name> <header value> and return the value.
template <class T>
T read_header(std::istream &stream, const std::string &name) {
	if (!stream.good()) {
		throw std::runtime_error("bad stream");
	}
	std::string header_name;
	T value;
	stream >> header_name >> value;
	// Check that the correct header and value were read
	if (!stream.good()) {
		throw std::runtime_error("not a header line");
	}
	if (header_name != name) {
		throw std::runtime_error(
				"invalid header name " + header_name + ". Expected " + name);
	}
	// Return the header value
	return value;
}

// Expect a certain character in a stream
struct expect {
    explicit expect(int c) : c(c) {}

    friend std::istream &operator>>(std::istream &stream, expect e) {
        using namespace std::literals;
        int x = stream.get();
        if (x != e.c) {
            throw std::runtime_error((std::ostringstream()
                                      << "expected char '" << e.c << "', received '" << x << "'"
                                      ).str());
        }
        return stream;
    }
private:
    int c;
};

Plant load_plant(const std::string &plant_path) {
	std::ifstream plant_istream(plant_path);
    if (!plant_istream.good()) {
        throw std::runtime_error("Unable to open " + plant_path);
    }

	Plant plant;

	// Read the number of vertices.
	plant.vertex_count = read_header<int>(plant_istream, "verts");
	// Allocate space to hold vertex data.
	plant.vertices = std::vector<Vertex>(plant.vertex_count);

	// Read all vertex data.
	for (int i = 0; i < plant.vertex_count; ++i) {
		Vertex vertex;
        plant_istream >> vertex.index >> expect(',')
                >> vertex.position[0] >> expect(',')
                >> vertex.position[1] >> expect(',')
                >> vertex.position[2] >> expect(',')
                >> vertex.radius >> expect(',')
                >> vertex.on_leaf;
        // Convert to z-up
        vertex.position = Matrix3d {
            { 1, 0, 0 },
            { 0, 0, 1 },
            { 0, -1, 0 },
        } * vertex.position;
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
	}
    plant.initDiffusion();
    plant.updateDiffusion(1.f);

	return plant;
}



