#pragma once

#include <string>
#include <QSettings>

struct IOConfig {
    std::string plant_file;
    std::string output_directory;
};

struct SimulationConfig {
    int frames;
};

struct WaterConfig {
    double dynamic_viscosity;
    double uniform_loss_rate;

    double time_between_frames;
    double integration_time_step;

    double elastic_steepness;
    double elastic_threshold;
};

struct RodsConfig {
    double density;
    double gravity;
    int iterations_per_frame;
    double iteration_time_step;
    int num_constraint_iterations;

    double min_youngs_modulus;
    double max_youngs_modulus;
    double torsion_modulus;
};

struct Config {
    IOConfig io;
    SimulationConfig simulation;
    WaterConfig water;
    RodsConfig rods;

    void init(const QSettings &settings);
};
