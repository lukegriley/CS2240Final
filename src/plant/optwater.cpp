#include <QCoreApplication>
#include <iostream>
#include <Eigen>

#include "Eigen/Sparse"

using namespace Eigen;

void printVector(const Vector3d& vec)
{
    std::cout << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")" << std::endl;
}

struct plant_node
{
    bool leaf;
    double volume;
    double loss_rate;
    double R;
    std::vector<int> neighbors;
};

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    std::cout << "\nhello world\n" << std::endl;



    // --- mesh representation ---

    // given plant geometry
    std::vector<Vector3f> vertices
    {
        Vector3f(0,0,0), // stem
        Vector3f(0,1,0), // stem
        Vector3f(0,2,0), // stem
        Vector3f(1,1,1)  // leaf
    };
    std::vector<Vector2i> edges
    {
        Vector2i(0, 1), // stem connection
        Vector2i(1, 2), // stem connection
        Vector2i(1, 3)  // leaf connection
    };

    // generate data structure for plant nodes from geometry
    std::vector<plant_node*> plant;
    for (int i = 0; i < vertices.size(); i++)
    {
        // object for node data structure
        plant_node* node = new plant_node;

        // leaf boolean... temp implementation
        node->leaf = vertices[i] == Vector3f(1,1,1);

        // find and store neighbors
        std::vector<int> neighbors;
        for (auto& edge : edges)
        {
            if (i != edge[0] && i != edge[1]) continue;
            neighbors.push_back(edge[0] != i? edge[0] : edge[1]);
        }
        node->neighbors = neighbors;

        // store object
        plant.push_back(node);
    }



    // --- water dispersion and evaporation initialization ---

    // variables - constant for entire plane
    double lv_main = 0.; // main structure loss rate
    double delta = 1.; // loss rate per surface area

    // per node variables
    for (auto& node: plant)
    {
        double rv = 1.; // radius of cylinder that approximates node
        double hv = 1.; // height of cylinder that approximates node
        node->volume = M_PI * rv * rv * hv; // volume

        double Av = rv * hv; // leaf surface area
        double lv_leaf = delta * Av; // leaf loss rate
        node->loss_rate = node->leaf? lv_leaf : lv_main;

        double mu = 1.; // dynamic viscosity
        double Rv = (M_PI * rv * rv * rv * rv) / (8 * mu * hv); // water flow resistance
        node->R = Rv;
    }

    // double thetav = 1.; // water amount contained in segment
    // double Pv = thetav / (rv * rv * M_PI * hv); // water pressure at the segment

    // D matrices
    SparseMatrix<double> Dv_inverse(vertices.size(), vertices.size()); // diagonal matrix with volumes
    Dv_inverse.reserve(1);
    for (int i = 0; i < vertices.size(); i++)
    {
        Dv_inverse.coeffRef(i,i) = 1 / plant[i]->volume; // volume at node...
    }
    Dv_inverse.makeCompressed();

    SparseMatrix<double> Dl(vertices.size(), vertices.size()); // diagonal matrix with loss rate
    Dl.reserve(1);
    for (int i = 0; i < vertices.size(); i++)
    {
        Dl.coeffRef(i,i) = plant[i]->loss_rate;
    }
    Dl.makeCompressed();

    // R matrix
    SparseMatrix<double> R(vertices.size(), vertices.size());
    for (int i = 0; i < vertices.size(); i++)
    {
        double sum = 0;
        for (auto neighbor : plant[i]->neighbors)
        {
            double resistance = (plant[i]->R + plant[neighbor]->R) / 2.;
            R.coeffRef(i, neighbor) = 1. / resistance;
            sum += 1. / resistance;
        }
        R.coeffRef(i,i) = -sum;
    }
    R.makeCompressed();

    // S combined matrix
    SparseMatrix<double> S = R * Dv_inverse - Dl;



    // --- continued initialization ---

    // inital value of water
    VectorXd theta(vertices.size());
    theta.setOnes();
    theta *= 50;

    // // section 4.1

    // spectral decomposition
    SelfAdjointEigenSolver<SparseMatrix<double>> spectral_decomp(S);
    assert(spectral_decomp.info() == Success);

    // stored value
    VectorXd theta_prime = spectral_decomp.eigenvectors().transpose() * theta;

    // --- iterated solver ---
    double dt = 0.01;
    double elapsed_time = 0.;
    while (elapsed_time <= 1.)
    {
        // fast computation
        VectorXd intermedidate_value = (spectral_decomp.eigenvalues() * elapsed_time).array().exp().matrix().asDiagonal() * theta_prime;
        VectorXd solution = spectral_decomp.eigenvectors() * intermedidate_value;

        std::cout << solution.transpose() << std::endl;

        // update total elapsed time
        elapsed_time += dt;
    }


    // section 4.2

    // variables
    double dt = 0.01;
    double elapsed_time = 0.;
    MatrixXd omega(vertices.size(), 6); omega.setOnes(); omega *= 50; // initial water values
    double alpha = 60. / 147.;
    VectorXd beta(6); beta << -10, 72, -225, 400, -450, 360; beta *= 1. / 147.;

    // precompute matrix decomposition
    SimplicialLLT<SparseMatrix<double>> solver(MatrixXd::Identity(vertices.size(), vertices.size()) - alpha * dt * S);

    // iterated solver
    while (elapsed_time <= 1.)
    {
        // b vector calculation
        VectorXd b = omega * beta;

        // fast computation
        VectorXd solution = solver.solve(b);

        std::cout << solution.transpose() << std::endl;

        // update omega
        omega.col(0) = omega.col(1);
        omega.col(1) = omega.col(2);
        omega.col(2) = omega.col(3);
        omega.col(3) = omega.col(4);
        omega.col(4) = omega.col(5);
        omega.col(5) = solution;

        // update total elapsed time
        elapsed_time += dt;
    }

    // --- export for python visualizations ---
    // todo...


    // return a.exec();
}
