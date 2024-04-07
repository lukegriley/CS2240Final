#include "Types.h"
#ifndef SYSTEM_H
#define SYSTEM_H
using namespace Eigen;
#include "collider.h"
#include "graphics/camera.h"

class System
{
public:
    System();
    void load(std::vector<Vector3d>& verts, std::vector<Vector4i>& tets, Matrix4d matrix);
    void loadModel(Matrix4d matrix) {
        m_model = matrix;
    }
    void loadCamera(Camera *_cam) {
        cam = _cam;
    }
    Matrix3d corotatedLinearElasticity();
    void zeroNodeForces();
    void computeParticleMasses();
    void update(double delta);
    void push(float negate);
    void addCollider(Collider col);
    void addColliders(std::vector<Collider> colliders);
    std::vector<Vector3d> eulerDerivEval();
    std::vector<Vector3d> getNodesPositions();
    bool br = false;
    void updateCameraPos(Vector3f o, Vector3f look) {
        cam_o=o;
        cam_look = look;
    }

private:
    void applyInternalForces();
    void applyCollisionForces();
    void checkCollisions();
    std::vector<Vector3d> rk4DerivEval(std::vector<Vector3d>& old_posAndVel, double delta);
    std::vector<Node> getNodesCopy();
    std::vector<Tet *> m_tets;
    std::vector<Node *> m_nodes;
    Vector3f cam_o;
    Vector3f cam_look;

    std::vector<Vector3d> addVectors(std::vector<Vector3d>& a, std::vector<Vector3d>& b) {
        assert(a.size()==b.size());
        std::vector<Vector3d> res(a.size());
        for(int i=0;i<a.size();i++) {
            res[i] = a[i] + b[i];
        }
    }

    std::vector<Vector3d> multiplyVectors(std::vector<Vector3d>& a, std::vector<Vector3d>& b) {
        assert(a.size()==b.size());
        std::vector<Vector3d> res(a.size());
        for(int i=0;i<a.size();i++) {
            res[i] = a[i] + b[i];
        }
    }

    Vector3d faceNormal(Tet *t, int node_i)
    {
        Node* n = t->nodes[node_i];
        std::vector<Node *> neighbors = {t->nodes[1],t->nodes[2],t->nodes[3]};
        if (node_i > 0) {
            neighbors[node_i - 1] = t->nodes[0];
        }

        Vector3d AB = neighbors[0]->pos-neighbors[1]->pos;
        Vector3d AC = neighbors[2]->pos-neighbors[1]->pos;

        Vector3d normal = AB.cross(AC).normalized();

        Vector3d AN = n->pos - neighbors[0]->pos;

        if (normal.dot(AN) >= 0) {
            return -normal;
        } else {
            return normal;
        }
    }
    void applyCollisionPenalty();
    void updateTetAreasAndNormals(Tet* t) {
        Vector3d side01(t->nodes[1]->pos - t->nodes[0]->pos);
        Vector3d side02(t->nodes[2]->pos - t->nodes[0]->pos);
        Vector3d side03(t->nodes[3]->pos - t->nodes[0]->pos);
        Vector3d side12(t->nodes[2]->pos - t->nodes[1]->pos);
        Vector3d side13(t->nodes[3]->pos - t->nodes[1]->pos);
        t->area012 = side01.cross(side02).norm() * 0.5;
        t->area013 = side01.cross(side03).norm() * 0.5;
        t->area023 = side02.cross(side03).norm() * 0.5;
        t->area123 = side12.cross(side13).norm() * 0.5;

        t->norm012 = faceNormal(t,3);
        t->norm013 = faceNormal(t,2);
        t->norm023 = faceNormal(t,1);
        t->norm123 = faceNormal(t,0);
    }

    Vector3d getWorldPos(Vector3d v) {
        Vector4d res = m_model * Vector4d(v[0],v[1],v[2],1);
        return Vector3d(res[0],res[1],res[2]);
    }

    Vector3d getMaterialPos(Vector3d v) {
        Vector4d res = m_model.inverse() * Vector4d(v[0],v[1],v[2],1);
        return Vector3d(res[0],res[1],res[2]);
    }

    Matrix4d m_model;

    //PARAMETERS
    double incompressibility =4e3;
    double rigidity = 4e3;
    double density = 12000.0;
    double phi = 1e3;
    double psi = 1e3;
    double gravity = -0.05;
    float normalFriction = 0.9f;
    float tangentialFriction = 0.2f;
    float pushForce = 1e3;
    std::vector<Vector3d> m_forcepushChanges;
    bool addForcePush = false;
    Camera *cam;


    std::vector<Collider> m_colliders;
};

#endif // SYSTEM_H
