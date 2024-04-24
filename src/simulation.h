#pragma once

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
private:
    Shape m_shape;
    Tree tree;
    RodRenderer renderer;

    Shape m_ground;
    void initGround();
};
