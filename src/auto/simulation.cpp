#include "simulation.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <plant/loader.h>

#include "export.h"

using namespace Eigen;

std::string filename(int n, int length) {
    std::stringstream sstr;
    sstr << "frame-"
         << std::setfill('0')
         << std::setw(length)
         << n
         << ".txt";
    return sstr.str();
}

Simulation::Simulation(const Config &config)
    : config(config) {

    // Initialize plant
    this->plant = plant::load(config.io.plant_file);

    this->water_plant.dynamic_viscosity = this->config.water.dynamic_viscosity;
    this->water_plant.delta = this->config.water.uniform_loss_rate;

    this->rod_tree.gravity = Vector3d {0, 0, -this->config.rods.gravity};
    this->rod_tree.num_iterations = this->config.rods.num_constraint_iterations;

}

void Simulation::run() {
    this->water_plant.initStructure(this->plant);
    this->water_plant.initDiffusion();
    this->rod_tree.init_from_plant(this->plant, config.rods.density);

    // Make output directory
    std::error_code ec;
    std::filesystem::path output_directory(this->config.io.output_directory);
    std::filesystem::create_directories(output_directory, ec);
    if (ec) {
        throw std::runtime_error("error creating output directory: " + ec.message());
    }


    for (int frame = 0; frame < config.simulation.frames; ++frame) {
        std::cout << "Running frame " << frame << std::endl;
        // Run the water simulation for some time
        for (double waterTime = 0;
             waterTime < config.water.time_between_frames;
             waterTime += config.water.integration_time_step) {
            this->water_plant.updateDiffusionDelta(config.water.integration_time_step);
        }
        std::cout << "Finished water" << std::endl;

        // Update the physics simulation
        for (int j = 0; j < this->water_plant.segments.size(); ++j) {
            const water::Segment &segment = this->water_plant.segments[j];
            double water_amount = segment.water_amt;
            double elasticity = 1 / (1 + std::exp(-this->config.water.elastic_steepness * (water_amount - this->config.water.elastic_threshold)));

            rod::Rod &rod = this->rod_tree.rods[j];
            double youngs_modulus = this->config.rods.min_youngs_modulus + (this->config.rods.max_youngs_modulus * elasticity);
            double torsion_modulus = this->config.rods.torsion_modulus;
            rod.set_material_parameters(youngs_modulus, torsion_modulus);
        }
        std::cout << "Updated parameters" << std::endl;

        // Run the physics simulation for some time
        for (int iteration = 0; iteration < this->config.rods.iterations_per_frame; ++iteration) {
            this->rod_tree.iterate(this->config.rods.iteration_time_step);
        }
        std::cout << "Finished physics" << std::endl;

        // Export the data
        std::cout << "Done with frame" << std::endl;

        std::filesystem::path out_path(this->config.io.output_directory);
        out_path /= filename(frame, 5);
        std::cout << "Writing to " << out_path << std::endl;
        std::ofstream out(out_path);
        if (!out.good()) {
            std::ostringstream msg;
            msg << "failed to open " << out_path << " for writing";
            throw std::runtime_error(msg.str());
        }
        out << this->rod_tree;
    }
}
