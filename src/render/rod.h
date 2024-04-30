#pragma once

#include <GL/glew.h>
#include <rods/rod.h>

#include "graphics/shader.h"
#include "graphics/shape.h"

// Renders rods as cubes
struct RodRenderer {
    RodRenderer();
    RodRenderer(const rod::Tree &tree);
    RodRenderer(const RodRenderer &renderer) = delete;

    void init(const rod::Tree &tree);
    void update(const rod::Tree &tree);
    void render(Shader *shader);
private:
    int rod_count;
    std::vector<Eigen::Affine3d> transformations;
    Shape m_cube;

    float uniform_scale = 3.f;
};


