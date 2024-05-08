#pragma once

#include <Eigen/Dense>

namespace rod {

// Apply interleaving order
int interleave(int n, int i);

// Use the given vector in quaternion math.
Eigen::Quaterniond as_quaternion(const Eigen::Vector3d &vec);

extern const Eigen::Matrix<double, 4, 3> as_quaternion_coeffs;

// Return a matrix M = M(v) such that Mx = cross(v, x).
Eigen::Matrix3d cross_product(const Eigen::Vector3d &v);

// Return the lower-right 3x3 of a the quaternion product.
Eigen::Matrix3d left_multiply_vec(const Eigen::Quaterniond &q);

// Return the lower-right 3x3 of a the quaternion product.
Eigen::Matrix3d left_multiply_vec_conjugate(const Eigen::Quaterniond &q);

// Represent left-multiplication by the given quaternion as a matrix.
Eigen::Matrix4d left_multiply(const Eigen::Quaterniond &q);

// Build G(q), according to equation 27.
// The differential of the quaternion q is given by G(q) * omega, where omega
// is the angular velocity.
Eigen::Matrix<double, 4, 3> angular_velocity_to_quaternion(const Eigen::Quaterniond &q);

}

