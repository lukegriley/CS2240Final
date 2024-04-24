#pragma once

#include <GL/glew.h>

#include "graphics/shader.h"
#include "graphics/shape.h"
#include "plant/plant.h"

// A PlantRenderer that draws cylinders for vessels
struct CylinderPlantRenderer {
    CylinderPlantRenderer();
    CylinderPlantRenderer(const Plant &plant);
    CylinderPlantRenderer(const CylinderPlantRenderer &renderer) = delete;

    void init(const Plant &plant);
    void update(const Plant &plant);
    void update_colors(const Plant &plant);
    void render(Shader *shader);
private:
    int node_count;
    std::vector<Eigen::Affine3d> transformations;
    std::vector<Eigen::Vector3f> colors;
    Shape m_cylinder;
    float uniform_scale = 10.f;
};


