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
    void render(Shader *shader);
private:
    int edge_count;
    std::vector<Eigen::Affine3d> transformations;
    Shape m_cylinder;
};


