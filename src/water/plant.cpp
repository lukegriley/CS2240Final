#include "water/plant.h"
#include <iostream>

using namespace std;
using namespace Eigen;

namespace water {

double Segment::height() const {
    return (this->tail_position - this->head_position).norm();
}

Segment createSegment(const plant::Plant &plant, const plant::Vertex &vertex) {
    return Segment {
        .index = vertex.index,
        .head_position = vertex.parent == -1 ? Vector3d::Zero() : plant.vertices[vertex.parent].tail_position,
        .tail_position = vertex.tail_position,
        .radius = vertex.radius,
        .on_leaf = vertex.on_leaf,
        .edges = vertex.edges,
    };
}

Edge createEdge(const plant::Edge edge) {
    return Edge {
        .segments = edge.vertices,
    };
}

void Plant::initStructure(const plant::Plant &basic_plant) {
    this->segments.resize(basic_plant.vertices.size());
    for (const plant::Vertex &vertex : basic_plant.vertices) {
        this->segments[vertex.index] = createSegment(basic_plant, vertex);
    }

    this->edges.resize(basic_plant.edges.size());
    for (const plant::Edge &edge : basic_plant.edges) {
        this->edges[edge.index] = createEdge(edge);
    }
}

void Plant::initDiffusion() {
    cout << "Starting diffusion init" <<endl;

    for (Segment &v : this->segments) {
        float height = v.height();

        // Compute loss rate as the approx. surface area covered by node * uniform loss rate (if on leaf)
        v.loss_rate = v.on_leaf ? height * v.radius * this->delta : 0.f;

        //Initialize theta (amount of water at a node) as its volume
        v.volume = M_PI * v.radius * v.radius * height;
        v.water_amt = v.volume;
    }

    //initialize D_V and D_l, diagonal matrices filled
    cout << "Initializing D_V and D_L" <<endl;
    SparseMatrix<double> sparse_DV_inverse(segments.size(), segments.size());
    SparseMatrix<double> sparse_DL(segments.size(), segments.size());
    sparse_DV_inverse.setZero(); sparse_DL.setZero();
    for(const Segment &v : this->segments) {
        sparse_DV_inverse.coeffRef(v.index, v.index) = 1. / v.volume;
        sparse_DL.coeffRef(v.index, v.index) = v.loss_rate;
    }
    sparse_DV_inverse.makeCompressed(); sparse_DL.makeCompressed();


    //compute pressure and resistance at each node
    cout << "Computing pressure and resistance at each node" <<endl;
    for(Segment &v : this->segments) {
        v.resistance = (M_PI * std::pow(v.radius,4))/(8.f*DYNAMIC_VISCOSITY*v.height());
    }

    //compute resistance at each segment (averaged between end nodes)
    for(Edge &e : this->edges) {
        e.resistance = (this->segments[e.segments[1]].resistance + this->segments[e.segments[0]].resistance)*0.5;
    }


    //initialize the symmetric matrix R
    SparseMatrix<double> R(this->segments.size(),this->segments.size());
    R.setZero();
    cout << "Initializing R" <<endl;

    for(Segment &v : this->segments) {
        double sum_resistance = 0.0;
        //go through each neighbor u
        for(int e_idx : v.edges) {
            const Edge &e = this->edges[e_idx];
            int u_idx = e.segments[0] == v.index ? e.segments[1] : e.segments[0];
            R.coeffRef(v.index, u_idx) = 1.0/e.resistance;
            R.coeffRef(u_idx, v.index) = 1.0/e.resistance;

            sum_resistance -= 1.0/e.resistance;
        }
        R.coeffRef(v.index,v.index) = sum_resistance;
    }

    m_S_sparse = R * sparse_DV_inverse - sparse_DL;

    // variables
    // initial water values!!!!!!!
    m_omega.resize(segments.size(), 6);
    for (const Segment &v : segments) {
        for (int col = 0; col < m_omega.cols(); ++col) {
            m_omega(v.index, col) = v.water_amt;
        }
    }
    m_alpha = 60. / 147.;
    m_beta.resize(6); m_beta << -10, 72, -225, 400, -450, 360; m_beta *= 1. / 147.;

    // cache matrix for fast computation
    // SparseMatrix<double> m_S_sparse = S.sparseView(); m_S_sparse.makeCompressed();
}

void Plant::updateDiffusionDelta(float dt)
{
    // sparse solver
    SparseMatrix<double> L(segments.size(), segments.size());
    L.setIdentity();
    L -= m_alpha * dt * m_S_sparse;
    L.makeCompressed();
    SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(L);
    solver.factorize(L);

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

    // Update water content in each node.
    for (Segment &segment : segments) {
        segment.water_amt = solution[segment.index];
    }
}

}
