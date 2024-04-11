#pragma once

#include <GL/glew.h>

#include "plant/plant.h"
#include "graphics/shader.h"

// A PlantRenderer draws lines
struct LinePlantRenderer {
    LinePlantRenderer() = default;
    LinePlantRenderer(const Plant &plant);
    LinePlantRenderer(const LinePlantRenderer &renderer) = delete;

	void init(const Plant &plant);
    void update(const Plant &plant);

    void render(Shader *shader) const;
private:
	int m_num_edges;

	GLuint m_vbo;
	GLuint m_ibo;
	GLuint m_vao;

    float uniform_scale = 10.f;
};


