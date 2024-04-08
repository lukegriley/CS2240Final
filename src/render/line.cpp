#include "line.h"

#include <iostream>

#include "util/unsupportedeigenthing/OpenGLSupport"

using namespace Eigen;

LinePlantRenderer::LinePlantRenderer() {
}

LinePlantRenderer::LinePlantRenderer(const Plant &plant) {
	this->init(plant);
}

void LinePlantRenderer::init(const Plant &plant) {
	this->m_num_edges = plant.edges.size();

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

	// Create a data vector with the edges
	std::vector<std::array<int, 2>> edges(plant.edges.size());
	for (int i = 0; i < plant.edges.size(); ++i) {
		edges[i] = plant.edges[i].vertices;
	}
	// Fill an IBO with edges
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
			sizeof(edges[0]) * plant.edges.size(),
			static_cast<const void *>(edges.data()),
			GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Set vertices
    this->update(plant);
}

void LinePlantRenderer::update(const Plant &plant) {
	// Create a data vector with the vertices
    std::vector<Vector3d> vertices(plant.vertices.size());
	for (int i = 0; i < plant.vertices.size(); ++i) {
		vertices[i] = plant.vertices[i].position;
	}
	// Fill the VBO with vertex positions
	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
	glBufferData(GL_ARRAY_BUFFER,
			sizeof(vertices[0]) * vertices.size(),
			static_cast<const void *>(vertices.data()),
			GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LinePlantRenderer::render() const {
	glBindVertexArray(m_vao);
	glDrawElements(
			GL_LINES,
            2 * m_num_edges,
			GL_UNSIGNED_INT,
			nullptr);
	glBindVertexArray(0);
}

