#include "cylinder.h"

#include <iostream>

#include "graphics/meshloader.h"
#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;

CylinderPlantRenderer::CylinderPlantRenderer() {
}

CylinderPlantRenderer::CylinderPlantRenderer(const Plant &plant) {
    this->init(plant);
}

void CylinderPlantRenderer::init(const Plant &plant) {
    // Load mesh topology
    std::vector<Vector3d> vertices;
    std::vector<Vector3i> faces;
    if (!MeshLoader::loadTriMesh("./meshes/stem.obj", vertices, faces)) {
        throw std::runtime_error("Unable to load cylinder mesh");
    }
    this->m_cylinder.init(vertices, faces);
    // Initialize each stem segment.
    this->update(plant);
}

void CylinderPlantRenderer::update(const Plant &plant) {
    this->edge_count = plant.edge_count;
    this->transformations.resize(plant.edge_count);
    // Set the transformation matrix as a uniform variable
    for (const Edge &edge : plant.edges) {
        const Vertex *vertices[2] {
            &plant.vertices[edge.vertices[0]],
            &plant.vertices[edge.vertices[1]],
        };
        double radius = vertices[1]->radius;
        Vector3d growth = vertices[1]->position - vertices[0]->position;

        // Construct a scaling matrix for the radius and length
        DiagonalMatrix<double, 3> scale(radius, radius, growth.norm());

        // Construct a rotation matrix that rotates +z to growth
        Vector3d rb[3];
        rb[2] = growth.normalized();
        rb[0] = rb[2].unitOrthogonal();
        rb[1] = rb[2].cross(rb[0]);
        Matrix3d rotation {
            { rb[0][0], rb[1][0], rb[2][0] },
            { rb[0][1], rb[1][1], rb[2][1] },
            { rb[0][2], rb[1][2], rb[2][2] },
        };
        assert((rotation * rotation.transpose() - Matrix3d::Identity()).cwiseAbs().maxCoeff() < 1e-3);

        // Construct a translation matrix to vertices[0]
        Translation3d translation(vertices[0]->position);

        this->transformations[edge.index] = translation * rotation * scale;
    }
}

void CylinderPlantRenderer::render(Shader *shader) {
    // For each edge, draw a cylinder.
    for (int i = 0; i < this->edge_count; ++i) {
        m_cylinder.setModelMatrix(this->transformations[i].cast<float>());
        // TODO: set the color as a uniform variable
        m_cylinder.draw(shader);
    }
}
