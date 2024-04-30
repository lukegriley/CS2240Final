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

//    std::vector<Vector3d> positions;
//    std::vector<double> masses;
//    std::vector<std::pair<int, int>> rods;
//    std::vector<double> radii;
//    Vector3d root { 0.1, 0, 0.2 };
//    double mass = 1e-6;

//    positions.push_back(root);
//    masses.push_back(mass);
//    int n = 100;
//    for (int i = 1; i < n; ++i) {
//        double angle = 0.4 * i;
//        positions.push_back(root + 0.05 * Vector3d {std::cos(angle) - 1, std::sin(angle), -0.05 * angle});
//        masses.push_back(mass);
//        rods.push_back(std::make_pair(i - 1, i));
//        radii.push_back(0.001);
//    }
//    // masses[masses.size() - 1] = mass * 1000;
//    tree.init_particles(positions, masses);
//    tree.init_orientations(rods, radii);

//    // Fix a rod
//    tree.rods[0].fixed = true;
//    tree.particles[0].fixed = true;
//    tree.particles[1].fixed = true;

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
    renderer.render(shader);
    m_ground.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

