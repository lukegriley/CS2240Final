#include "cylinder.h"

#include <iostream>

#include "graphics/meshloader.h"
#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;
using namespace water;

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
    this->node_count = plant.segments.size();
    this->transformations.resize(this->node_count);
    // Set the transformation matrix as a uniform variable
    for (const Segment &segment : plant.segments) {
        double radius = segment.radius;
        Vector3d growth = segment.tail_position - segment.head_position;

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
        Translation3d translation(segment.head_position);

        this->transformations[segment.index] = translation * rotation * scale;
    }

    this->update_colors(plant);
}

void CylinderPlantRenderer::update_colors(const Plant &plant) {
    // Initialize colors
    this->colors.resize(plant.segments.size());
    for (const Segment &segment : plant.segments) {
        float density = segment.water_amt / segment.volume;
        // Clamp density
        density = std::min(std::max(density, 0.0f), 1.99f);
        // Split 0-2 into 4 pieces for red-green and green-blue.
        float H = 2 * density;
        switch (static_cast<int>(std::floor(H))) {
        case 0:
            this->colors[segment.index] = Vector3f(1.0f, H - 0, 0.0f);
            break;
        case 1:
            this->colors[segment.index] = Vector3f(1.0f - (H - 1), 1.0f, 0.0f);
            break;
        case 2:
            this->colors[segment.index] = Vector3f(0.0f, 1.0f, H - 2);
            break;
        case 3:
            this->colors[segment.index] = Vector3f(0.0f, 1.0f - (H - 3), 1.0f);
            break;
        }
    }
}

void CylinderPlantRenderer::render(Shader *shader) {
    // Convert to y-up
    const Matrix3f toYUp {
        { 1, 0, 0 },
        { 0, 0, 1 },
        { 0, -1, 0 },
    };

    // For each edge, draw a cylinder.
    for (int i = 0; i < this->node_count; ++i) {
        Matrix3f scaleConstant = Matrix3f::Identity() * this->uniform_scale;
        m_cylinder.setModelMatrix(scaleConstant * toYUp * this->transformations[i].cast<float>());
        shader->setUniform("color", this->colors[i]);
        m_cylinder.draw(shader);
    }
}
