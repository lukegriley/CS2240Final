#pragma once

#include <GL/glew.h>

#include "water/plant.h"
#include "graphics/shader.h"

// A PlantRenderer draws lines
struct LinePlantRenderer {
    LinePlantRenderer() = default;
    LinePlantRenderer(const water::Plant &plant);
    LinePlantRenderer(const LinePlantRenderer &renderer) = delete;

    void init(const water::Plant &plant);
    void update(const water::Plant &plant);

    void render(Shader *shader) const;
private:
    int m_num_nodes;

	GLuint m_vbo;
	GLuint m_ibo;
	GLuint m_vao;

    float uniform_scale = 10.f;
};


