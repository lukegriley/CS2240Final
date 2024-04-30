#pragma once

#include <GL/glew.h>

#include "graphics/shader.h"
#include "graphics/shape.h"
#include "water/plant.h"

// A PlantRenderer that draws cylinders for vessels
struct CylinderPlantRenderer {
    CylinderPlantRenderer();
    CylinderPlantRenderer(const water::Plant &plant);
    CylinderPlantRenderer(const CylinderPlantRenderer &renderer) = delete;

    void init(const water::Plant &plant);
    void update(const water::Plant &plant);
    void update_colors(const water::Plant &plant);
    void render(Shader *shader);
private:
    int node_count;
    std::vector<Eigen::Affine3d> transformations;
    std::vector<Eigen::Vector3f> colors;
    Shape m_cylinder;
    float uniform_scale = 10.f;
};


