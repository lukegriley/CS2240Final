#include "plant/loader.h"

#include <fstream>

using namespace Eigen;

namespace plant {


// Read a line formatted as <header name> <header value> and return the value.
template <class T>
static T read_header(std::istream &stream, const std::string &name) {
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

Plant load(const std::string &path) {
    std::ifstream plant_istream(path);
    if (!plant_istream.good()) {
        throw std::runtime_error("Unable to open " + path);
    }

    Plant plant;

    // Read the number of vertices.
    int vertex_count = read_header<int>(plant_istream, "verts");
    // Allocate space to hold vertex data.
    plant.vertices = std::vector<Vertex>(vertex_count);

    // Read all vertex data.
    for (int i = 0; i < vertex_count; ++i) {
        Vertex vertex;
        plant_istream >> vertex.index >> expect(',')
                >> vertex.tail_position[0] >> expect(',')
                >> vertex.tail_position[1] >> expect(',')
                >> vertex.tail_position[2] >> expect(',')
                >> vertex.radius >> expect(',')
                >> vertex.on_leaf;
        // Read optional arguments
        if (plant_istream.get() == ',') {
            plant_istream >> vertex.fixed;
        }

        assert(vertex.index == i);
        plant.vertices[i] = vertex;
    }

    // Read the number of edges
    int edge_count = read_header<int>(plant_istream, "edges");
    // Allocate space for all edges
    plant.edges = std::vector<Edge>(edge_count);
    // Read all edge data.
    for (int i = 0; i < edge_count; ++i) {
        Edge edge { .index = i };
        plant_istream >> edge.vertices[0] >> expect(',')
                >> edge.vertices[1];

        // Add the edge to the plant
        plant.edges[i] = edge;
        // Add the edge to all vertices
        for (int j = 0; j < 2; ++j) {
            plant.vertices[edge.vertices[j]].edges.push_back(i);
        }
        // Record the parent
        if (plant.vertices[edge.vertices[1]].parent != -1) {
            throw std::runtime_error("vertex has two parents");
        } else {
            plant.vertices[edge.vertices[1]].parent = edge.vertices[0];
        }
        // Record the child
        plant.vertices[edge.vertices[0]].children.push_back(edge.vertices[1]);
    }

    return plant;
}

}

