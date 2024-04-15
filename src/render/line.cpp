#include "line.h"

#include <iostream>

#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;

LinePlantRenderer::LinePlantRenderer(const Plant &plant) {
    this->init(plant);
}

void LinePlantRenderer::init(const Plant &plant) {
    this->m_num_nodes = plant.vertices.size();

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

void LinePlantRenderer::update(const Plant &plant) {
    // Create a data vector with the vertices
    std::vector<Vector3f> vertices(2 * plant.vertices.size());
    for (int i = 0; i < plant.vertices.size(); ++i) {
        vertices[2 * i] = plant.vertices[i].head_position.cast<float>();
        vertices[2 * i + 1] = plant.vertices[i].tail_position.cast<float>();
    }
    // Fill the VBO with vertex positions
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(vertices[0]) * vertices.size(),
            static_cast<const void *>(vertices.data()),
            GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LinePlantRenderer::render(Shader *shader) const {
    Matrix4f model = Transform<float, 3, Affine>(Scaling(uniform_scale)).matrix();
    shader->setUniform("model", model);
    Eigen::Matrix3f inverseTransposeModel = 1.f / uniform_scale * Matrix3f::Identity();
    shader->setUniform("inverseTransposeModel", inverseTransposeModel);

    glBindVertexArray(m_vao);
    glDrawArrays(GL_LINES,
                 0,
                 2 * m_num_nodes);
    glBindVertexArray(0);
}

