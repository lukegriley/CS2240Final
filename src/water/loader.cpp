#include "water/loader.h"

#include <fstream>

#include "plant/loader.h"

using namespace Eigen;

namespace water {

void Loader::load_S_decomp(const std::string &plant_path, Plant &plant){
    std::ifstream plant_istream(plant_path);
    if (!plant_istream.good()) {
        throw std::runtime_error("Unable to open " + plant_path);
    }
    int phi_size = Loader::read_header<int>(plant_istream, "phi");
    plant.phi = Eigen::MatrixXf::Zero(phi_size,phi_size);

    for (int i = 0; i < phi_size; ++i) {
        for(int j=0;j<phi_size;j++) {
            float coeff;
            plant_istream >> coeff >> expect(' ');
            plant.phi(i,j) = coeff;
        }
        expect('\n');
    }

    // int lambda_size = Loader::read_header<int>(plant_istream, "lambda");
    // plant.lambda = Eigen::MatrixXf::Zero(lambda_size,lambda_size);

    // for (int i = 0; i < lambda_size; ++i) {
    //     float coeff;
    //     plant_istream >> coeff >> expect('\n');
    //     plant.lambda(i,i) = coeff;
    // }

}

}
