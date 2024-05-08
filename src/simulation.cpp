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
    m_plant.dynamic_viscosity = this->dynamic_viscosity;
    m_plant.delta = this->uniform_loss_rate;
    m_plant.initStructure(plant);
    m_plant.initDiffusion();
    m_plantRenderer.init(m_plant);

    // Create the rod system
    double density = this->density;
    tree.init_from_plant(plant, density);
//    for (rod::direct::Constraint &constraint : tree.constraints) {
//        constraint.set_material_parameters(1e2, 1e2);
//    }
//     tree.rods[2].position += Vector3d {-1, 0, 0};

    tree.gravity = this->gravity;
    renderer.init(tree);
}

void Simulation::update(double seconds)
{
    m_plant.updateDiffusionDelta(seconds);
    this->m_plantRenderer.update_colors(m_plant);

//    for (rod::direct::Constraint &constraint : tree.constraints) {
//        int j = constraint.rods[1];
//        const water::Segment &segment = m_plant.segments[j];
//        double water_amount = segment.water_amt; // segment.volume;
//        double elasticity = 1 / (1 + std::exp(-elastic_steepness * (water_amount - elastic_threshold)));

//        double youngs_modulus = 1e8 + (1e10 * elasticity);
//        double torsion_modulus = 79000000000;
//        double stiffness = 1; //1e-4; //1e11; // 1e-4;
//        // constraint.set_material_parameters(stiffness * youngs_modulus, stiffness * torsion_modulus);
//    }

    for (rod::Rod &rod : tree.rods) {
        int j = rod.index;
        const water::Segment &segment = m_plant.segments[j];
        double water_amount = segment.water_amt;
        double elasticity = 1 / (1 + std::exp(-elastic_steepness * (water_amount - elastic_threshold)));

        double youngs_modulus = this->min_youngs_modulus + (this->max_youngs_modulus * elasticity);
        double torsion_modulus = this->torsion_modulus;
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

