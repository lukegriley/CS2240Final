#include "rodutils.h"

#include <cassert>

using namespace Eigen;

namespace rod {

int interleave_even(int n, int i) {
    assert(n % 2 == 0);
    if (i % 2 == 0) {
        return i;
    } else {
        return n - i;
    }
}

// Apply interleaving order
int interleave(int n, int i) {
    if (n % 2 == 0) {
        return interleave_even(n, i);
    } else {
        if (i == 0) {
            return 0;
        } else {
            return 1 + interleave_even(n - 1, i - 1);
        }
    }
}

// Use the given vector in quaternion math.
Quaterniond as_quaternion(const Vector3d &vec) {
    return Quaterniond(0, vec[0], vec[1], vec[2]);
}

// Use the given vector in quaternion math.
const Matrix<double, 4, 3> as_quaternion_coeffs {
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 0, 0, 1 },
    { 0, 0, 0 },
};


// Return a matrix M = M(v) such that Mx = cross(v, x).
Matrix3d cross_product(const Vector3d &v) {
    return Matrix3d {
        { 0, -v[2], v[1] },
        { v[2], 0, -v[0] },
        { -v[1], v[0], 0 },
    };
}

// Return the lower-right 3x3 of a the quaternion product.
Matrix3d left_multiply_vec(const Quaterniond &q) {
    return q.w() * Matrix3d::Identity() + cross_product(q.vec());
}

// Return the lower-right 3x3 of a the quaternion product.
Matrix3d left_multiply_vec_conjugate(const Quaterniond &q) {
    return q.w() * Matrix3d::Identity() - cross_product(q.vec());
}

Matrix4d left_multiply(const Quaterniond &q) {
    Matrix4d qm;
    // Place the real coordinate at the end, like Eigen.
    qm(3, 3) = q.w();
    qm.block<1, 3>(3, 0) = -q.vec().transpose();
    qm.block<3, 1>(0, 3) = q.vec();
    qm.block<3, 3>(0, 0) = left_multiply_vec(q);
    return qm;
}

// Build G(q), according to equation 27.
// The differential of the quaternion q is given by G(q) * omega, where omega
// is the angular velocity.
Matrix<double, 4, 3> angular_velocity_to_quaternion(const Quaterniond &q) {
    return 0.5 * left_multiply(q) * as_quaternion_coeffs;

    Matrix<double, 4, 3> G;
    // https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/master/PositionBasedDynamics/PositionBasedElasticRods.cpp
    // Lines 658 to 667.
    // The paper says to place this in row 0, but we place it in row 3 because
    // Eigen stores the scalar at the end.
    G.row(3) = -q.vec().transpose();
    // Note that this is distinct from left_multiply_vec.
    // left_multiply_vec is in the paper.
    // left_multiply_vec_conjugate is in the code.
    G.block<3, 3>(0, 0) = left_multiply_vec(q);
    return 0.5 * G;
}

}
