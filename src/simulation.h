#pragma once

#include "plant/plant.h"
#include "render/cylinder.h"
#include "graphics/shape.h"
#include "render/rod.h"

class Shader;

class Simulation
{
public:
    Simulation();

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    void toggleWire();

    double elastic_threshold = 3e-8;
    double elastic_steepness = 1e8;
private:
    Shape m_shape;
    water::Plant m_plant;
    CylinderPlantRenderer m_plantRenderer;
    rod::Tree tree;
    RodRenderer renderer;

    Shape m_ground;
    void initGround();
};
