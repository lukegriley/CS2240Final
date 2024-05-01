#include "simulation.h"

#include <iostream>

#include <plant/loader.h>
#include <water/loader.h>

#include "graphics/meshloader.h"

using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init()
{

    std::string plant_file = "./data/plants/plant000.txt";
    plant::Plant plant = plant::load(plant_file);

    // Initialize plant and renderer
    m_plant = water::Loader::load_plant(plant_file);
    m_plantRenderer.init(m_plant);

    // Create the rod system
    double density = 7800;
    tree.init_from_plant(plant, density);

    tree.gravity = Vector3d {0, 0, -1};
    renderer.init(tree);

    // Set the number of timesteps
    tree.num_bend_twist_steps = 1;
}

void Simulation::update(double seconds)
{
    m_plant.updateDiffusionDelta(seconds);
    this->m_plantRenderer.update_colors(m_plant);

    tree.iterate(seconds);
    renderer.update(tree);
}

void Simulation::draw(Shader *shader)
{
    this->m_plantRenderer.render(shader);
    renderer.render(shader);
    m_ground.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

