#include <Eigen/Dense>

#ifndef COLLIDER_H
#define COLLIDER_H
#include "Types.h"
using namespace Eigen;


enum ColliderType {
    COLL_PLANE,
    COLL_SPHERE,
    COLL_AAB,
};



using namespace Eigen;

class Collider
{
public:
    Collider(ColliderType _type);
    IntInfo checkIntersection(Vector3d point);
    void loadPlane(int _axis, float _height, float _min1, float _max1, float _min2, float _max2);
    void loadSphere(Vector3d center, float _radius);
    void loadAAB(float _min1, float _max1, float _min2, float _max2, float _min3, float _max3);
    ColliderType m_type;
private:
    IntInfo planeIntersection(Vector3d point);
    IntInfo sphereIntersection(Vector3d point);
    IntInfo AABIntersection(Vector3d point);
    int plane_axis;
    int plane_height;
    float min1;//for AAB and plane
    float max1;//for AAB and plane
    float min2;//for AAB and plane
    float max2;//for AAB and plane
    float min3;//for AAB
    float max3;//for AAB
    Vector3d center;//for spheres
    float radius;//for spheres
    float epsilon = .1;
};

#endif // COLLIDER_H
