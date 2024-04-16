#pragma once

#include "plant.h"

class Writer {
public:

    static void write_plant(const std::string &plant_path, const Plant &plant);
    static void write_S_decomp(const std::string &S_decomp_path, const Plant &plant);
};
