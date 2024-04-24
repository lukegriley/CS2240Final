#pragma once

#include <GL/glew.h>

#include "graphics/shader.h"
#include "graphics/shape.h"
#include "rods/rod.h"

// Renders rods as cubes
struct RodRenderer {
    RodRenderer();
    RodRenderer(const Tree &tree);
    RodRenderer(const RodRenderer &renderer) = delete;

    void init(const Tree &tree);
    void update(const Tree &tree);
    void render(Shader *shader);
private:
    int rod_count;
    std::vector<Eigen::Affine3d> transformations;
    Shape m_cube;
};


