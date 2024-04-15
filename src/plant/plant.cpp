#include "plant.h"
#include <iostream>

using namespace std;
using namespace Eigen;

/*
 * NOTE from Luke 4/10:
 * This has been implemented according to the analytical solution (5) in section 3.2
 * As the paper mentions in section 4, this is a naive method to compute theta_t
 * since the matrix multiplication is very expensive.
 *
 * so TODO:we need to implement section 4.1 and 4.2 to make this more efficient
 */

void Plant::initDiffusion(bool precompute) {
    cout << "Starting diffusion init" <<endl;
    this->theta_0 = Eigen::MatrixXf(this->vertices.size(),1);

    for(Edge e : this->edges) {
        Vertex *v0 = &this->vertices[e.vertices[0]];
        Vertex *v1 = &this->vertices[e.vertices[1]];
        float height = (v1->position - v0->position).norm();

        // Compute loss rate as the approx. surface area covered by node * uniform loss rate (if on leaf)
        v0->loss_rate = v0->on_leaf ? height * v0->radius * this->delta : 0.f;

        //Initialize theta (amount of water at a node) as its volume
        v0->volume = M_PI * v0->radius * v0->radius * height;
        v0->water_amt = v0->volume;
        this->theta_0.coeffRef(v0->index,0) = v0->volume;
    }


    // //initialize D_V and D_l, diagonal matrices filled
    // cout << "Initializing D_V and D_L" <<endl;
    // this->D_V = Eigen::MatrixXf(this->vertices.size(),this->vertices.size());
    // this->D_V.setZero();
    // this->D_l = Eigen::MatrixXf(this->vertices.size(),this->vertices.size());
    // this->D_l.setZero();
    // for(Vertex &v : this->vertices) {
    //     this->D_V.coeffRef(v.index,v.index) = v.volume;
    //     this->D_l.coeffRef(v.index,v.index) = v.loss_rate;
    // }

    //initialize D_V and D_l, diagonal matrices filled
    cout << "Initializing D_V and D_L" <<endl;
    SparseMatrix<double> sparse_DV_inverse(vertices.size(), vertices.size());
    SparseMatrix<double> sparse_DL(vertices.size(), vertices.size());
    sparse_DV_inverse.setZero(); sparse_DL.setZero();
    for(Vertex &v : this->vertices) {
        sparse_DV_inverse.coeffRef(v.index,v.index) = 1. / v.volume;
        sparse_DL.coeffRef(v.index,v.index) = v.loss_rate;
    }
    sparse_DV_inverse.makeCompressed(); sparse_DL.makeCompressed();


    //compute pressure and resistance at each node
    cout << "Computing pressure and resistance at each node" <<endl;
    for(Edge e : this->edges) {
        Vertex *v0 = &this->vertices[e.vertices[0]];
        Vertex *v1 = &this->vertices[e.vertices[1]];
        v0->resistance = (M_PI * std::pow(v0->radius,4))/(8.f*DYNAMIC_VISCOSITY*(v1->position - v0->position).norm());
    }

    //compute resistance at each segment (averaged between end nodes)
    for(Edge &e : this->edges) {
        e.resistance = (this->vertices[e.vertices[1]].resistance + this->vertices[e.vertices[0]].resistance)*0.5;
    }


    //initialize the symmetric matrix R
    // Eigen::MatrixXf R(this->vertices.size(),this->vertices.size());
    SparseMatrix<double> R(this->vertices.size(),this->vertices.size());
    R.setZero();
    cout << "Initializing R" <<endl;

    for(Vertex &v : this->vertices) {
        float sum_resistance = 0.0;
        //go through each neighbor u
        for(int e_idx : v.edges) {
            Edge &e = this->edges[e_idx];
            int u_idx = e.vertices[0]==v.index ? e.vertices[1] : e.vertices[0];
            R.coeffRef(v.index,u_idx) = 1.0/e.resistance;
            R.coeffRef(u_idx,v.index) = 1.0/e.resistance;

            sum_resistance -= e.resistance;
        }
        R.coeffRef(v.index,v.index) = sum_resistance;
    }

    // cout << "Computing S" <<endl;//this step will take a minute or two
    // S = R * this->D_V.inverse() - this->D_l;
    m_S_sparse = R * sparse_DV_inverse - sparse_DL;

    // // For testing, skip the slow precomputations.
    // if (!precompute)
    //     return;

    // //a start on 4.1
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(S);
    // if (eigensolver.info() != Eigen::Success) {
    //     cerr << "Spectral decomp of S matrix failed" <<endl;
    // } else {
    //     Eigen::VectorXf eigenvalues = eigensolver.eigenvalues();

    //     this->phi = eigensolver.eigenvectors();
    //     this->lambda.diagonal() = eigenvalues;
    // }


    // 4.2 precomputation for fast solver

    // variables
    MatrixXd m_omega(vertices.size(), 6); m_omega.setOnes(); m_omega *= 50; // initial water values!!!!!!!
    double m_alpha = 60. / 147.;
    VectorXd m_beta(6); m_beta << -10, 72, -225, 400, -450, 360; m_beta *= 1. / 147.;

    // cache matrix for fast computation
    // SparseMatrix<double> m_S_sparse = S.sparseView(); m_S_sparse.makeCompressed();


}



void Plant::updateDiffusion(float time) {
    // cout << "Computing theta at t="<<time <<endl;
    // Eigen::MatrixXf theta_t = (this->S*time).array().exp().matrix() * this->theta_0;
    // for(int i=0;i<this->vertices.size();i++) {
    //     this->vertices[i].water_amt = theta_t.coeffRef(i,0);
    // }
    // cout << "Update finished" <<endl;

    cout << "computation moved to updateDiffusionDelta!!" << std::endl;
}

void Plant::updateDiffusionDelta(const float dt)
{
    // sparse solver
    SimplicialLLT<SparseMatrix<double>> solver(MatrixXd::Identity(vertices.size(), vertices.size()) - m_alpha * dt * m_S_sparse);

    // b vector calculation
    VectorXd b = m_omega * m_beta;

    // fast computation
    VectorXd solution = solver.solve(b);

    // update omega
    m_omega.col(0) = m_omega.col(1);
    m_omega.col(1) = m_omega.col(2);
    m_omega.col(2) = m_omega.col(3);
    m_omega.col(3) = m_omega.col(4);
    m_omega.col(4) = m_omega.col(5);
    m_omega.col(5) = solution;
}

