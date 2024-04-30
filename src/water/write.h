#pragma once

#include "plant.h"

namespace water {

class Writer {
public:

    static void write_S_decomp(const std::string &S_decomp_path, const Plant &plant);
};

}
