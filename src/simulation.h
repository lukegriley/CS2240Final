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

    double dynamic_viscosity = 9e-7;
    double uniform_loss_rate = 1e3;

    double density = 7800;
    Eigen::Vector3d gravity {0, 0, -1e-4};
    double min_youngs_modulus = 1;
    double max_youngs_modulus = 1e10;
    double elastic_threshold = 3e-7;
    double elastic_steepness = 1e8;
    double torsion_modulus = 79000000000;
private:
    Shape m_shape;
    water::Plant m_plant;
    CylinderPlantRenderer m_plantRenderer;
    rod::Tree tree;
    RodRenderer renderer;
//    rod::direct::Tree tree;
//    DirectRodRenderer renderer;

    Shape m_ground;
    void initGround();
};
