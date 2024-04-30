#include "water/write.h"

#include <fstream>

namespace water {

void Writer::write_S_decomp(const std::string &S_decomp_path, const Plant &plant) {
    std::ofstream S_ostream(S_decomp_path);
    if (!S_ostream.good()) {
        throw std::runtime_error("Unable to open " + S_decomp_path);
    }

    S_ostream << "phi " << plant.phi.rows() <<" " << plant.phi.cols() << std::endl;
    for (int i = 0; i < plant.phi.rows(); ++i) {
        for (int j = 0; j < plant.phi.cols(); ++j) {
            S_ostream << plant.phi(i, j) << " ";
        }
        S_ostream << std::endl;
    }
    S_ostream << "lambda " << plant.lambda.rows() << std::endl;
    for (int i = 0; i < plant.lambda.rows(); ++i) {
        S_ostream << plant.phi(i, i);
        S_ostream << std::endl;
    }

    S_ostream.close();
}

}
