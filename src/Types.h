#ifndef TYPES_H
#define TYPES_H
#include <Eigen/Dense>
struct Tet;
struct Node;
using namespace Eigen;



struct IntInfo {
    bool hit = false;
    Vector3d surfacePoint = Vector3d(0,0,0);
    Vector3d surfaceNormal = Vector3d(0,0,0);
    float penetration;
    Vector3d force;
};

struct Tet {
    Node *nodes[4];
    double volume;
    Matrix3d beta;
    Matrix<double,3,4> x;
    double area012;
    double area013;
    double area123;
    double area023;
    Vector3d norm012;
    Vector3d norm013;
    Vector3d norm123;
    Vector3d norm023;
    Matrix<double, 3,4> internalForceChange;
    Matrix<double, 3, 12> N;
};

struct Node {
    double mass;
    Vector3d pos;
    Vector3d force;
    Vector3d velocity;
    int index;
    IntInfo i;
};


namespace std {
template <>
struct hash<Vector3i> {
    size_t operator()(const Vector3i& v) const {
        size_t hashValue = hash<int>()(v[0]);
        hashValue ^= hash<int>()(v[1]) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
        hashValue ^= hash<int>()(v[2]) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
        return hashValue;
    }
};
}





#endif // TYPES_H
