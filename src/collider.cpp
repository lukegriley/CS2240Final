#include "collider.h"
#include <iostream>
using namespace Eigen;
using namespace std;


Collider::Collider(ColliderType _type) {
    m_type = _type;
}

void Collider::loadPlane(int _axis, float _height, float _min1, float _max1, float _min2, float _max2) {
    if (m_type!=ColliderType::COLL_PLANE) {
        cerr << "Collider is not a plane, cannot load geometry." <<std::endl;
    }
    plane_axis = 1;
    plane_height = _height;
    min1 = _min1;
    max1 = _max1;
    min2 = _min2;
    max2 = _max2;
}
void Collider::loadSphere(Vector3d _center, float _radius) {
    if (m_type!=ColliderType::COLL_SPHERE) {
        cerr << "Collider is not a sphere, cannot load geometry." <<std::endl;
    }
    center = _center;
    radius = _radius;
}

void Collider::loadAAB(float _min1, float _max1, float _min2, float _max2, float _min3, float _max3) {
    if (m_type!=ColliderType::COLL_AAB) {
        cerr << "Collider is not an axis-aligned box, cannot load geometry." <<std::endl;
    }
    min1 = _min1;
    max1 = _max1;
    min2 = _min2;
    max2 = _max2;
    min3 = _min3;
    max3 = _max3;
}

IntInfo Collider::checkIntersection(Vector3d point) {
    switch(m_type) {
    case COLL_AAB:
        return AABIntersection(point);
    case COLL_PLANE:
        return planeIntersection(point);
    case COLL_SPHERE:
        return sphereIntersection(point);
    }
}

IntInfo Collider::planeIntersection(Vector3d point) {
    // if(point[plane_axis] < plane_height + epsilon && point[plane_axis] > plane_height - epsilon) {

    // }
    IntInfo res;
    //for now, just treating like an infinite ground plane
    int axis1 = (plane_axis + 1 )%3;
    int axis2 = (plane_axis + 2)%3;
    if(point[plane_axis] < plane_height && point[axis1] > min1 && point[axis1] < max1
        && point[axis2] > min2 && point[axis2] < max2) {
        res.penetration = abs(plane_height - point[plane_axis]);
        point[plane_axis] = plane_height;
        res.surfacePoint = point;
        Vector3d normal(0,0,0);
        normal[plane_axis] = 1.f;
        res.surfaceNormal = normal;
        res.hit = true;

        return res;
    }
    return res;
}
IntInfo Collider::sphereIntersection(Vector3d point) {
    double dist = (point-center).norm();
    IntInfo res;
    if(dist < radius) {
        res.hit = true;
        Vector3d direction = (point - center).normalized();
        res.surfacePoint = center + radius * direction;
        res.surfaceNormal = direction;
        res.penetration = radius-dist;
    }
    return res;
}


IntInfo Collider::AABIntersection(Vector3d point) {

}
