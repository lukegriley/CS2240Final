#include "directrod.h"

#include <iostream>
#include <GL/glew.h>

#include "graphics/meshloader.h"
#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;
using namespace rod::direct;

DirectRodRenderer::DirectRodRenderer() {
}

DirectRodRenderer::DirectRodRenderer(const Tree &tree) {
    this->init(tree);
}

void DirectRodRenderer::init(const Tree &tree) {
    // Load mesh topology
    std::vector<Vector3d> vertices;
    std::vector<Vector3i> faces;
    if (!MeshLoader::loadTriMesh("./meshes/cube.obj", vertices, faces)) {
        throw std::runtime_error("Unable to load cube mesh");
    }
    this->m_cube.init(vertices, faces);
    // Initialize each stem segment.
    this->update(tree);
}

void DirectRodRenderer::update(const Tree &tree) {
    this->rod_count = tree.rods.size();
    this->transformations.resize(this->rod_count);
    // Set the transformation matrix as a uniform variable
    for (const Rod &rod : tree.rods) {
        double radius = rod.radius;
        Vector3d growth = rod.direction();

        // Construct a scaling matrix for the radius and length
        DiagonalMatrix<double, 3> scale(radius, radius, growth.norm());

        const auto &rotation = rod.orientation;

        // Construct a translation matrix to the rod center
        Translation3d translation(rod.position);

        this->transformations[rod.index] = translation * rotation * scale;
    }
}

void DirectRodRenderer::render(Shader *shader) {
    // Convert to y-up
    const Matrix3f toYUp {
        { 1, 0, 0 },
        { 0, 0, 1 },
        { 0, -1, 0 },
    };

    // For each edge, draw a cylinder.
    for (int i = 0; i < this->rod_count; ++i) {
        m_cube.setModelMatrix(toYUp * uniform_scale * this->transformations[i].cast<float>());
        // TODO: set the color as a uniform variable
        m_cube.draw(shader);
    }
}
