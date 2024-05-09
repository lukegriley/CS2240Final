#pragma once

#include <plant/plant.h>
#include <rods/rod.h>
#include <water/plant.h>

#include "config.h"

struct Simulation {
    Config config;
    plant::Plant plant;
    water::Plant water_plant;
    rod::Tree rod_tree;

    Simulation(const Config &config);
    void run();
    void export_water(int frame);
    void export_physics(int frame);
};
