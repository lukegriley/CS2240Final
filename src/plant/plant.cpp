#include "plant.h"
#include <iostream>
#include <fstream>
#include "loader.h"
#include "write.h"

using namespace std;


/*
 * NOTE from Luke 4/10:
 * This has been implemented according to the analytical solution (5) in section 3.2
 * As the paper mentions in section 4, this is a naive method to compute theta_t
 * since the matrix multiplication is very expensive.
 *
 * so TODO:we need to implement section 4.1 and 4.2 to make this more efficient
 */

void Plant::initDiffusion() {
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
        this->theta_0.coeffRef(v0->index,0) = v0->volume;
    }


    //initialize D_V and D_l, diagonal matrices filled
    cout << "Initializing D_V and D_L" <<endl;
    this->D_V = Eigen::MatrixXf(this->vertices.size(),this->vertices.size());
    this->D_V.setZero();
    this->D_l = Eigen::MatrixXf(this->vertices.size(),this->vertices.size());
    this->D_l.setZero();
    for(Vertex &v : this->vertices) {
        this->D_V.coeffRef(v.index,v.index) = v.volume;
        this->D_l.coeffRef(v.index,v.index) = v.loss_rate;
    }


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

    if(!load_S) {
        //initialize the symmetric matrix R
        Eigen::MatrixXf R(this->vertices.size(),this->vertices.size());
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
        cout << "Computing D V inv" <<endl;//this step will take a minute or two
        auto D_V_inv = this->D_V.inverse();
        cout << "Computing S" <<endl;//this step will take a minute or two
        S = R * D_V_inv - this->D_l;

        this->phi = Eigen::MatrixXf(this->vertex_count,this->vertex_count);
        this->lambda = Eigen::MatrixXf(this->vertex_count,this->vertex_count);
        cout << "Computing S decomp" <<endl;//this step will take a minute or two

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(S);
        if (eigensolver.info() != Eigen::Success) {
            cerr << "Spectral decomp of S matrix failed" <<endl;
        } else {
            Eigen::VectorXf eigenvalues = eigensolver.eigenvalues();

            this->phi = eigensolver.eigenvectors();
            this->lambda.diagonal() = eigenvalues;
        }

        Writer::write_S_decomp("./data/S_decomp/plant0.txt",*this);
    } else {
        Loader::load_S_decomp("./data/S_decomp/plant0.txt",*this);
    }

    // cout << "Computing zeta_0" <<endl;
    // this->zeta_0 = this->phi.colPivHouseholderQr().solve(this->theta_0);
    cout<<"init finished" <<endl;

}



void Plant::updateDiffusion(float time) {
    cout << "Computing theta at t="<<time <<endl;
    // Eigen::MatrixXf theta_t = (this->S*time).array().exp().matrix() * this->theta_0;
    Eigen::MatrixXf theta_t = this->phi * ((this->lambda*time).array().exp().matrix() * this->zeta_0);

    for(int i=0;i<this->vertices.size();i++) {
        this->vertices[i].water_amt = theta_t.coeffRef(i,0);
    }
    cout << "Update finished" <<endl;
}






