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

    tree.gravity = Vector3d {0, 0, -10};
    renderer.init(tree);

    // Set the number of timesteps
    tree.num_bend_twist_steps = 1;
}

void Simulation::update(double seconds)
{
    m_plant.delta = 3;
    m_plant.updateDiffusionDelta(10 * seconds);
    this->m_plantRenderer.update_colors(m_plant);

    for (int i = 0; i < m_plant.segments.size(); ++i) {
        const water::Segment &segment = m_plant.segments[i];
        rod::Rod &rod = tree.rods[i];
        double water_amount = segment.water_amt; // segment.volume;
        double elasticity = 1 / (1 + std::exp(-elastic_steepness * (water_amount - elastic_threshold)));

        double youngs_modulus = 1 + (1e10 * elasticity);
        double torsion_modulus = 79000000000;
        rod.set_material_parameters(youngs_modulus, torsion_modulus);
    }

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

