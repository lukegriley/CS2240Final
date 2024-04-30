#include "line.h"

#include <iostream>

#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;

LinePlantRenderer::LinePlantRenderer(const water::Plant &plant) {
    this->init(plant);
}

void LinePlantRenderer::init(const water::Plant &plant) {
    this->m_num_nodes = plant.segments.size();

    // Doesn't always work. We will develop a better renderer soon
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(3);

    // Generate a VBO, an IBO, and a VAO
    glGenBuffers(1, &m_vbo);
    glGenBuffers(1, &m_ibo);
    glGenVertexArrays(1, &m_vao);

    // Bind the VBO and IBO to the VAO
    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Set vertices
    this->update(plant);
}

void LinePlantRenderer::update(const water::Plant &plant) {
    // Create a data vector with the vertices
    std::vector<Vector3f> segments(2 * plant.segments.size());
    for (int i = 0; i < plant.segments.size(); ++i) {
        segments[2 * i] = plant.segments[i].head_position.cast<float>();
        segments[2 * i + 1] = plant.segments[i].tail_position.cast<float>();
    }
    // Fill the VBO with vertex positions
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(segments[0]) * segments.size(),
            static_cast<const void *>(segments.data()),
            GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LinePlantRenderer::render(Shader *shader) const {
    // Convert to y-up
    const Matrix3f toYUp {
        { 1, 0, 0 },
        { 0, 0, 1 },
        { 0, -1, 0 },
    };

    Matrix4f model = Matrix4f::Identity();
    model.topLeftCorner(3, 3) = toYUp * Scaling(uniform_scale);
    shader->setUniform("model", model);
    Eigen::Matrix3f inverseTransposeModel = 1.f / uniform_scale * Matrix3f::Identity();
    shader->setUniform("inverseTransposeModel", inverseTransposeModel);

    glBindVertexArray(m_vao);
    glDrawArrays(GL_LINES,
                 0,
                 2 * m_num_nodes);
    glBindVertexArray(0);
}

