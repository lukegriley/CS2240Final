#include "simulation.h"

#include <plant/loader.h>

Simulation::Simulation(const Config &config)
    : config(config) {

    // Initialize plant
    this->plant = plant::load(config.io.plant_file);
    this->water_plant.initStructure(this->plant);
    this->rod_tree.init_from_plant(this->plant, config.rods.density);
}

void Simulation::run() {
    for (int frame = 0; frame < config.simulation.frames; ++frame) {
        // Run the water simulation for some time
        for (double waterTime = 0;
             waterTime < config.water.time_between_frames;
             waterTime += config.water.integration_time_step) {
            this->water_plant.updateDiffusionDelta(config.water.integration_time_step);
        }
        // Initialize the physics simulation
        for (int i = 0; i < this->water_plant.segments.size(); ++i) {
            const water::Segment &segment = this->water_plant.segments[i];
            rod::Rod &rod = this->rod_tree.rods[i];

            double density = segment.water_amt / segment.volume;
            double elasticity = 1 / (1 + exp(-config.water.elastic_steepness * (density - config.water.elastic_threshold)));
            rod.set_material_parameters(elasticity, config.rods.torsion_modulus);
        }
        // Run the physics simulation for some time
    }
}
