#include "waterball.h"

#include <iostream>

#include "graphics/meshloader.h"
#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;

WaterBallPlantRenderer::WaterBallPlantRenderer(const Plant &plant) {
    this->init(plant);
}

void WaterBallPlantRenderer::init(const Plant &plant) {
    // Load mesh topology
    std::vector<Vector3d> vertices;
    std::vector<Vector3i> faces;
    if (!MeshLoader::loadTriMesh("./meshes/sphere.obj", vertices, faces)) {
        throw std::runtime_error("Unable to load cylinder mesh");
    }
    this->m_ball.init(vertices, faces);
    // Initialize each stem segment.
    this->update(plant);
}

void WaterBallPlantRenderer::update(const Plant &plant) {
    this->vertex_count = plant.vertex_count;
    this->transformations.resize(plant.vertex_count);
    // Set the transformation matrix as a uniform variable
    for (const Vertex &vertex : plant.vertices) {
        // Construct a scaling matrix for the water volume
        float radius = std::pow(vertex.volume, 0.33f);
        DiagonalMatrix<double, 3> scale(radius, radius, radius);
        // Construct a translation matrix to vertices[0]
        Translation3d translation(vertex.position);
        this->transformations[vertex.index] = translation * scale;
    }
    // Set colors
    this->colors.resize(plant.vertex_count);
    for (const Vertex &vertex : plant.vertices) {
        // Construct a scaling matrix for the water volume
        float density = vertex.water_amt / vertex.volume;
        // Clamp density
        density = std::min(std::max(density, 0.0f), 1.99f);
        // Split 0-2 into 4 pieces for red-green and green-blue.
        float H = 2 * density;
        switch (static_cast<int>(std::floor(H))) {
        case 0:
            this->colors[vertex.index] = Vector3f(1.0f, H - 0, 0.0f);
            break;
        case 1:
            this->colors[vertex.index] = Vector3f(1.0f - (H - 1), 1.0f, 0.0f);
            break;
        case 2:
            this->colors[vertex.index] = Vector3f(0.0f, 1.0f, H - 2);
            break;
        case 3:
            this->colors[vertex.index] = Vector3f(0.0f, 1.0f - (H - 3), 1.0f);
            break;
        }
    }
}

void WaterBallPlantRenderer::render(Shader *shader) {
    // For each vertex, draw a ball.
    for (int i = 0; i < this->vertex_count; ++i) {
        Matrix3f scaleConstant = Matrix3f::Identity() * this->uniform_scale;
        m_ball.setModelMatrix(scaleConstant * this->transformations[i].cast<float>());
        shader->setUniform("color", this->colors[i]);
        m_ball.draw(shader);
    }
}
