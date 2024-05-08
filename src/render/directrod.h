#pragma once

#include <GL/glew.h>
#include <rods/direct.h>

#include "graphics/shader.h"
#include "graphics/shape.h"

// Renders rods as cubes
struct DirectRodRenderer {
    DirectRodRenderer();
    DirectRodRenderer(const rod::direct::Tree &tree);
    DirectRodRenderer(const DirectRodRenderer &renderer) = delete;

    void init(const rod::direct::Tree &tree);
    void update(const rod::direct::Tree &tree);
    void render(Shader *shader);
private:
    int rod_count;
    std::vector<Eigen::Affine3d> transformations;
    Shape m_cube;

    float uniform_scale = 10.f;
};


