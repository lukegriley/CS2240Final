#include "system.h"
#include <iostream>

System::System() {

}

using namespace std;

void System::load(std::vector<Vector3d>& verts, std::vector<Vector4i>& tets, Matrix4d matrix) {
    m_nodes = std::vector<Node *>(verts.size());
    m_tets = std::vector<Tet *>(tets.size());
    m_forcepushChanges = std::vector<Vector3d>(verts.size());
    m_model = matrix;
    for(int i = 0;i<verts.size();i++) {
        m_nodes[i] = new Node{.pos = verts[i], .index = i};
    }

    for(int i = 0;i<tets.size();i++) {
        Tet *t = new Tet {.nodes = {
                m_nodes[tets[i][0]],
                m_nodes[tets[i][1]],
                m_nodes[tets[i][2]],
                m_nodes[tets[i][3]]}};
        m_tets[i] = t;

        Matrix3d beta;
        beta.col(0) = (t->nodes[0]->pos-t->nodes[3]->pos);
        beta.col(1) = (t->nodes[1]->pos-t->nodes[3]->pos);
        beta.col(2) = (t->nodes[2]->pos-t->nodes[3]->pos);
        // beta.col(0) = (getWorldPos(t->nodes[0]->pos) - getWorldPos(t->nodes[3]->pos));
        // beta.col(1) = (getWorldPos(t->nodes[1]->pos) - getWorldPos(t->nodes[3]->pos));
        // beta.col(2) = (getWorldPos(t->nodes[2]->pos) - getWorldPos(t->nodes[3]->pos));
        t->beta = beta.inverse();

        //shape function
        Matrix<double, 3, 12> N;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 3; ++j) {
                N(j, i * 3 + j) = 1.0 - t->nodes[i]->pos[j];
                N(j, (i + 1) * 3 + j) = t->nodes[i]->pos[j];
            }
            t->x.col(i) = Vector3d(t->nodes[i]->pos[0],t->nodes[i]->pos[1],t->nodes[i]->pos[2]);
        }
        t->N = N;


        // t->nodes[0]->pos += Vector3d(0,0.5,0);
        updateTetAreasAndNormals(t);
    }
    computeParticleMasses();
};



void System::addCollider(Collider col) {
    m_colliders.push_back(col);
}

void System::zeroNodeForces() {
    for(Node *n : m_nodes) {
        n->force = Vector3d(0.0,0.0,0.0);
    }
}

void System::computeParticleMasses() {
    for(Node *n : m_nodes) {
        n->mass = 0.0;
        n->velocity = Vector3d(0.0,0.0,0.0);
    }
    for(Tet *t : m_tets) {
        Eigen::Vector3d a = t->nodes[0]->pos;
        Eigen::Vector3d b = t->nodes[1]->pos;
        Eigen::Vector3d c = t->nodes[2]->pos;
        Eigen::Vector3d d = t->nodes[3]->pos;

        double volume = std::abs((1.0 / 6.0) * (a - d).dot((b - d).cross(c - d)));
        double nodeMass = volume * density * 0.25;
        t->volume = volume;

        t->nodes[0]->mass += nodeMass;
        t->nodes[1]->mass += nodeMass;
        t->nodes[2]->mass += nodeMass;
        t->nodes[3]->mass += nodeMass;
    }
}



void System::update(double delta) {
    if(br) {
        int x = 0;
    }

    // auto maxNode = std::max_element(m_nodes.begin(), m_nodes.end(), [](Node* a, Node* b) {
    //     return a->force.norm() < b->force.norm();
    // });
    // if(maxNode !=m_nodes.end()) {
    //     cout << (*maxNode)->force.norm() <<endl;
    // }

    std::vector<Vector3d> old_posAndVel = std::vector<Vector3d>(m_nodes.size() * 2);
    for(int i = 0; i < m_nodes.size(); i++) {
        old_posAndVel[2 * i] = m_nodes[i]->pos;
        old_posAndVel[2 * i + 1] = m_nodes[i]->velocity;
    }

    std::vector<Vector3d> k1 = eulerDerivEval();
    for(int i = 0; i < m_nodes.size(); i++) {
        m_nodes[i]->pos = old_posAndVel[2*i] + k1[2*i] * delta;
        m_nodes[i]->velocity = old_posAndVel[2*i + 1] + k1[2*i + 1] * delta;
    }
    std::vector<Vector3d> k2 = eulerDerivEval();
    for(int i = 0; i < m_nodes.size(); i++) {
        m_nodes[i]->pos = old_posAndVel[2*i] + k2[2*i] * delta;
        m_nodes[i]->velocity = old_posAndVel[2*i + 1] + k2[2*i + 1] * delta;
    }
    std::vector<Vector3d> k3 = eulerDerivEval();
    for(int i = 0; i < m_nodes.size(); i++) {
        m_nodes[i]->pos = old_posAndVel[2*i] + k3[2*i] * delta;
        m_nodes[i]->velocity = old_posAndVel[2*i + 1] + k3[2*i + 1] * delta;
    }


    std::vector<Vector3d> k4 = eulerDerivEval();
    for(int i = 0; i < m_nodes.size(); i++) {   //delta/6 * k1 * 2k2 * 2k3 * k4 + old
        m_nodes[i]->pos = old_posAndVel[2*i] + (k1[i*2] + 2*(k2[i*2] + k3[i*2]) + k4[i*2]) * (delta/6.0);
        m_nodes[i]->velocity = old_posAndVel[2*i + 1] + (k1[i*2 + 1] + 2*(k2[i*2 + 1] + k3[i*2 + 1]) + k4[i*2 + 1]) * (delta/6.0);
    }

    applyCollisionForces();
    if(addForcePush) {
        addForcePush = false;
    }
}

std::vector<Vector3d> System::eulerDerivEval() {
    // zero out forces + apply gravity
#pragma omp parallel for
    for(Node *n : m_nodes) {
        n->force = Vector3d(0.0,gravity,0.0) * n->mass;
        // n->force = Vector3d(0.0,0.0,0.0);
    }

    // applyCollisionPenalty();
    if(addForcePush) {
        for(int i = 0; i < m_nodes.size(); i++) {
            m_nodes[i]->force += m_forcepushChanges[i];
        }
    }

    //internal forces: strain and stress
    applyInternalForces();
    for(int i=0;i<m_tets.size();i++){
        for(int k=0;k<4;k++) {
            m_tets[i]->nodes[k]->force += Vector3d(m_tets[i]->internalForceChange.col(k));
        }
    }

    std::vector<Vector3d> destination;
    for(Node *n : m_nodes) {
        destination.push_back(n->velocity);
        if(n->mass<=0) {
            cerr << "node mass <= 0" << endl;
        }
        destination.push_back(n->force/n->mass);
    }
    return destination;
}

void System::applyCollisionPenalty() {
    for(Node *n : m_nodes) {
        n->i = m_colliders[0].checkIntersection(getWorldPos(n->pos));
        Vector3d collForce = 6e3 * n->i.surfaceNormal;
        n->force += collForce;
    }
}

void System::applyInternalForces() {
    // #pragma omp parallel for default(none) shared(m_tets) private(i, t, P, F, VP, V, strain, strain_rate, elastic_stress, viscous_damping, stress, force012, force013, force023, force123, force0, force1, force2, force3)
#pragma omp for parallel
    for(int i=0;i<m_tets.size();i++){
        Tet *t = m_tets[i];
        // updateTetAreasAndNormals(t);
        Matrix3d P;
        P.col(0) = (getWorldPos(t->nodes[0]->pos) - getWorldPos(t->nodes[3]->pos));
        P.col(1) = (getWorldPos(t->nodes[1]->pos) - getWorldPos(t->nodes[3]->pos));
        P.col(2) = (getWorldPos(t->nodes[2]->pos) - getWorldPos(t->nodes[3]->pos));
        // P.col(0) = (t->nodes[0]->pos-t->nodes[3]->pos);
        // P.col(1) = (t->nodes[1]->pos-t->nodes[3]->pos);
        // P.col(2) = (t->nodes[2]->pos-t->nodes[3]->pos);
        Matrix3d F = P * t->beta;
        if(F.isIdentity()) {
            F = Matrix3d::Identity();
        }

        Matrix3d VP;
        VP.col(0) = (t->nodes[0]->velocity - t->nodes[3]->velocity);
        VP.col(1) = (t->nodes[1]->velocity - t->nodes[3]->velocity);
        VP.col(2) = (t->nodes[2]->velocity - t->nodes[3]->velocity);
        Matrix3d V = VP * t->beta;
        if(V.isIdentity()) {
            V = Matrix3d::Identity();
        }

        Matrix3d strain = F.transpose() * F - Matrix3d::Identity();
        Matrix3d strain_rate = F.transpose() * V + V.transpose() * F;

        Matrix3d elastic_stress = incompressibility * Matrix3d::Identity() * strain.trace()
                                  + 2.0 * rigidity * strain;

        Matrix3d viscous_damping = phi * Matrix3d::Identity() * strain_rate.trace() +
                                   2.0 * psi * strain_rate;

        Matrix3d stress = elastic_stress + viscous_damping;

        Vector3d force012 = F*stress * t->area012 * t->norm012;
        Vector3d force013 = F*stress * t->area013 * t->norm013;
        Vector3d force023 = F*stress * t->area023 * t->norm023;
        Vector3d force123 = F*stress * t->area123 * t->norm123;

        // t->nodes[0]->force += (-1.0/3.0)*(force012 + force013 + force023);
        // t->nodes[1]->force += (-1.0/3.0)*(force012 + force013 + force123);
        // t->nodes[2]->force += (-1.0/3.0)*(force012 + force123 + force023);
        // t->nodes[3]->force += (-1.0/3.0)*(force123 + force013 + force023);
        Vector3d force0 = (-1.0/3.0)*(force012 + force013 + force023);
        Vector3d force1 = (-1.0/3.0)*(force012 + force013 + force123);
        Vector3d force2 = (-1.0/3.0)*(force012 + force123 + force023);
        Vector3d force3 = (-1.0/3.0)*(force123 + force013 + force023);
        t->internalForceChange.col(0) = force0;
        t->internalForceChange.col(1) = force1;
        t->internalForceChange.col(2) = force2;
        t->internalForceChange.col(3) = force3;

        if(force012.norm() > 100000000) {
            // cout << "over" <<endl;
        }

    }
}

void System::applyCollisionForces() {
#pragma omp parallel for
    for(Node *n : m_nodes) {
        for(Collider c: m_colliders) {
            IntInfo i = c.checkIntersection(getWorldPos(n->pos));
            if(i.hit) {
                n->pos = getMaterialPos(i.surfacePoint) + Vector3d(0,0.01,0);
                Vector3d normalComp = i.surfaceNormal * n->velocity.dot(i.surfaceNormal);
                Vector3d tangentialComp = n->velocity - normalComp;


                n->velocity = -normalComp * normalFriction + tangentialComp * tangentialFriction;

            }
        }
    }
}

Matrix3d System::corotatedLinearElasticity() {
    for (Tet* t : m_tets) {
        // Compute the strain-displacement matrix B for the tetrahedron

        // Compute linearized strain tensor
        // Assuming small deformations, so the linearized strain tensor is epsilon = 0.5 * (grad(u) + (grad(u))^T)

        // Compute linearized stress tensor
        // Assuming linear elastic material, so stress = stiffnessMatrix * strain

        // Compute the tangent stiffness matrix
        // Assuming the corotated linear elasticity method, so K = B^T * C * B, where C is the material stiffness tensor

        // Accumulate the contribution to the global stiffness matrix
        // stiffnessMatrix += t->beta.transpose() * t->volume; // Assuming beta is the material stiffness tensor
    }


}

void System::push(float negate) {
    //camera origin o
    //camera ray dir d
    //Vi vertex pos
    //(Vi - o) • d * forceScalar * d * attentuation
    Vector3f o = cam_o;
    Vector3f d = cam_look;
    for(int i=0;i<m_nodes.size();i++) {
        Vector3f v = getWorldPos(m_nodes[i]->pos).cast<float>();
        float distance = (v-o).norm();
        float att = 1.f/(.4f + 0.1f*distance + 0.1f*distance * distance);
        float dotted = std::clamp((v-o).dot(d)-(float)cos(1.f*M_PI/3.f),0.f,1.f);
        float force_scalar = pushForce * (1.f-dotted);
        m_forcepushChanges[i] = negate * (double)force_scalar * (v-o).normalized().cast<double>();
        // m_forcepushChanges[i] = (double)pushForce * (v-o).normalized().cast<double>();
    }
    cout << "push" <<endl;
    addForcePush = true;
}

std::vector<Vector3d> System::getNodesPositions() {
    std::vector<Vector3d> res(m_nodes.size());
    for(int i=0;i<m_nodes.size();i++) {
        res[i] = m_nodes[i]->pos;
    }
    return res;
}
