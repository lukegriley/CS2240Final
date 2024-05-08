#include "direct.h"

#include <iostream>
#include <stack>

#include <omp.h>
#include <plant/plant.h>

#include "rodutils.h"

using namespace Eigen;

namespace rod::direct {



/////// Math utils


Rod::Rod(int index,
         const Vector3d &position,
         double radius,
         double length,
         double density,
         bool fixed)
    : index(index),
      position(position),
      radius(radius),
      length(length),
      density(density),
      fixed(fixed),
      parent(-1),
      orientation(0, 0, 0, 0),
      velocity(0, 0, 0),
      angular_velocity(0, 0, 0)
{
}

double Rod::volume() const {
    return std::numbers::pi * this->radius * this->radius * this->length;
}

double Rod::mass() const {
    return this->density * this->volume();
}

Vector3d Rod::direction() const  {
    return this->length * (this->orientation * Vector3d(0, 0, 1));
}

double Rod::area_inertia() const {
    return 0.25 * std::numbers::pi * radius * radius * radius * radius;
}

Matrix<double, 6, 6> Rod::mass_matrix(const Quaterniond &orientation) const {
    Matrix<double, 6, 6> mass = Matrix<double, 6, 6>::Identity();

    // If the rod is fixed, pretend to have a very large mass.
    if (this->fixed) {
        mass *= this->mass() * 1e5;
        return mass;
    }

    // Set the mass in the first 3 entries
    mass.block<3, 3>(0, 0) *= this->mass();
    // Set the inertia in the last 3 entries
    double area_inertia = this->area_inertia();
    double inertia = this->density * this->length * area_inertia;
//     inertia = 1.0 / 12 * this->mass() * this->length * this->length;
    double roll_inertia = 2 * inertia;
    //    double roll_inertia = 0.5 * this->mass() * this->radius * this->radius;

    DiagonalMatrix<double, 3> local_inertia(inertia, inertia, roll_inertia);
    // Convert the inertia from local coordinates to world coordinates
    // https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/master/PositionBasedDynamics/PositionBasedElasticRods.cpp
    // Lines 708 to 712.
    Matrix3d rotation = orientation.toRotationMatrix();
    Matrix3d world_inertia = rotation * local_inertia * rotation.transpose();
    mass.block<3, 3>(3, 3) = world_inertia;
    return mass;
}



//// Darboux vector


Darboux::Darboux(const Quaterniond &q1, const Quaterniond &q2, double length)
    : q {q1, q2},
      length(length)
{
}

Vector3d Darboux::vec() const {
    return 2 / this->length * (this->q[0].conjugate() * this->q[1]).vec();
}

// Equations 10 and 11
Matrix<double, 3, 4> Darboux::differentiate(int rod_index) const {
    int other_index = 1 - rod_index;
    Matrix<double, 3, 4> d = as_quaternion_coeffs.transpose() * left_multiply(this->q[other_index].conjugate());

    int sign = (rod_index == 1) ? 1 : -1;
    return sign * 2 / this->length * d;
}


/// Constraints

Constraint::Constraint(int index, const Rod &rod1, const Rod &rod2)
    : index(index),
      rods {rod1.index, rod2.index},
      length {rod1.length, rod2.length},
      radii {rod1.radius, rod2.radius},
      initial_darboux(Darboux(rod1.orientation, rod2.orientation, 0.5 * (rod1.length + rod2.length)).vec()),
      compliance {1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8}
{
}

double Constraint::median_length() const {
    return 0.5 * (this->length[0] + this->length[1]);
}

void Constraint::set_material_parameters(double youngs_modulus, double torsion_modulus) {
    double stretch_stiffness = 1e9;
    double r = 0.5 * (this->radii[0] + this->radii[1]);
    double area_inertia = 0.25 * std::numbers::pi * r * r * r * r;

    // Refer to line 1120 of
    // https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/master/PositionBasedDynamics/PositionBasedElasticRods.cpp
    // for the multiplication by the length.
    double bend_stiffness =  youngs_modulus * area_inertia;
    double twist_stiffness = torsion_modulus * area_inertia * 2;
    DiagonalMatrix<double, 6> material = {
        stretch_stiffness,
        stretch_stiffness,
        stretch_stiffness,
        bend_stiffness,
        bend_stiffness,
        twist_stiffness,
    };
    this->compliance = material.inverse();
}

DiagonalMatrix<double, 6> Constraint::compliance_step(double dt) const {
    return DiagonalMatrix<double, 6>(this->compliance.diagonal() / (dt * dt));
}


Vector<double, 6> Constraint::evaluate(
        const std::vector<Vector3d> &positions,
        const std::vector<Quaterniond> &orientations) const {
    const Vector3d p1 = 0.5 * this->length[0] * Vector3d(0, 0, 1);
    const Vector3d p2 = 0.5 * this->length[1] * Vector3d(0, 0, -1);
    const Vector3d &x1 = positions[this->rods[0]];
    const Vector3d &x2 = positions[this->rods[1]];
    const Quaterniond &q1 = orientations[this->rods[0]];
    const Quaterniond &q2 = orientations[this->rods[1]];

    Vector3d new_darboux = Darboux(q1, q2, this->median_length()).vec();
    Vector3d initial_darboux = this->initial_darboux;

    Vector<double, 6> constraint;
    constraint.head(3) = q1 * p1 + x1 - q2 * p2 - x2;
    constraint.tail(3) = new_darboux - initial_darboux; //initial_darboux;
    return constraint;
}

Matrix<double, 6, 6> Constraint::jacobian(
        int rod_index,
        const Vector3d &p,
        const Quaterniond &q,
        const Darboux &darboux) const {
    Matrix<double, 6, 6> jacobian;
    Matrix<double, 4, 3> G = angular_velocity_to_quaternion(q);

    // Negate some terms for the child
    int sign = (rod_index == 0) ? 1 : -1;

    // Set the top-left 3x3 to the identity
    jacobian.block<3, 3>(0, 0) = sign * Matrix3d::Identity();
    assert((jacobian.block<3, 3>(0, 0).allFinite()));
    // Set the top-right to the cross product of the constraint point.
    jacobian.block<3, 3>(0, 3) = -sign * cross_product(q * p);

    // https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
    Quaterniond p_quat = as_quaternion(p);
    Quaterniond i(0, 1, 0, 0);
    Quaterniond j(0, 0, 1, 0);
    Quaterniond k(0, 0, 0, 1);
    Matrix4d left;
    left.col(0) = (p_quat * q * i).conjugate().coeffs() - (p_quat * q * i).coeffs();
    left.col(1) = (p_quat * q * j).conjugate().coeffs() - (p_quat * q * j).coeffs();
    left.col(2) = (p_quat * q * k).conjugate().coeffs() - (p_quat * q * k).coeffs();
    left.col(3) = (p_quat * q).coeffs() - (p_quat * q).conjugate().coeffs();
    jacobian.block<3, 3>(0, 3) = sign * as_quaternion_coeffs.transpose() * left * G;
    assert((jacobian.block<3, 3>(0, 3).allFinite()));
    // Set the bottom left to be zero
    jacobian.block<3, 3>(3, 0).setZero();
    assert((jacobian.block<3, 3>(3, 0).allFinite()));
    // Set the bottom right to the derivative of the darboux vector.
    jacobian.block<3, 3>(3, 3) = darboux.differentiate(rod_index) * G;
    assert((jacobian.block<3, 3>(3, 3).allFinite()));

    return jacobian;
}



void Tree::init(
        const std::vector<Vector3d> &positions,
        const std::vector<double> &radii,
        const std::vector<Vector3d> &directions,
        const std::vector<double> &densities,
        const std::vector<bool> &fixed,
        const std::vector<std::pair<int, int>> &edges) {

    int m = positions.size();
    assert(m == radii.size());
    assert(m == directions.size());
    assert(m == densities.size());
    assert(m == fixed.size());
    this->rods.resize(m);

    // Compute particles and radius for each rod
    for (int j = 0; j < m; ++j) {
        this->rods[j] = Rod(j,
                            positions.at(j),
                            radii.at(j),
                            directions.at(j).norm(),
                            densities.at(j),
                            fixed.at(j));

        assert(this->rods[j].length > 0);
    }

    // Compute quaternions for all rods
    std::vector<int> parents;
    parents.resize(m, -1);
    for (const std::pair<int, int> &edge : edges) {
        assert(parents.at(edge.second) == -1);
        parents.at(edge.second) = edge.first;
    }
    this->compute_rest_orientation(directions, parents);

    // Initialize constraints from plant edges
    int n = edges.size();
    this->constraints.resize(n);
    for (int i = 0; i < n; ++i) {
        const std::pair<int, int> &edge  = edges[i];
        this->constraints[i] = Constraint(
                    i,
                    this->rods[std::get<0>(edge)],
                    this->rods[std::get<1>(edge)]);
    }

    // Initialize H
    const int rod_size = 6;
    const int constraint_size = 6;
    int size = rod_size * m + constraint_size * n;
    this->H.resize(size, size);
}

void Tree::compute_rest_orientation(const std::vector<Vector3d> &directions,
                                    const std::vector<int> &parents) {
    int m = this->rods.size();
    std::vector<bool> visited(m, false);
    std::vector<Quaterniond> orientations(m);
    std::stack<int> to_visit;

    for (int j = 0; j < m; ++j) {
        to_visit.push(j);
    }
    while (to_visit.size() > 0) {
        int next = to_visit.top();
        if (visited[next]) {
            to_visit.pop();
            continue;
        }
        int parent_index = parents[next];

        // Set root, if the rod has no parent
        if (parent_index == -1) {
            orientations[next].setFromTwoVectors(Vector3d{0, 0, 1}, directions[next]);
            orientations[next].normalize();
            to_visit.pop();
            visited[next] = true;
            continue;
        }

        // Visit parent first
        if (!visited[parent_index]) {
            to_visit.push(parent_index);
            continue;
        }

        // Initialize quaternion
        orientations[next] = Quaterniond::FromTwoVectors(directions[parent_index], directions[next]) * orientations[parent_index];
        orientations[next].normalize();
        to_visit.pop();
        visited[next] = true;
    }
    for (int j = 0; j < m; ++j) {
        this->rods[j].orientation = orientations[j];
    }
}

void Tree::init_from_plant(const plant::Plant &plant, double density) {
    using std::numbers::pi;

    // Initialize first position to be zero

    std::vector<Vector3d> positions(plant.vertices.size());
    std::vector<double> radii(plant.vertices.size());
    std::vector<Vector3d> directions(plant.vertices.size());
    std::vector<double> densities(plant.vertices.size());
    std::vector<bool> fixed(plant.vertices.size());
    std::vector<std::pair<int, int>> edges(plant.edges.size());

    // Define rods from vertices
    for (const plant::Vertex &vertex : plant.vertices) {
        positions.at(vertex.index) = 0.5 * (vertex.head_position(plant) + vertex.tail_position);
        radii.at(vertex.index) = vertex.radius;
        directions.at(vertex.index) = vertex.direction(plant);
        densities.at(vertex.index) = density;
        fixed.at(vertex.index) = vertex.fixed || vertex.parent == -1;
    }

    // Define constraints from edges
    for (const plant::Edge &edge : plant.edges) {
        edges.at(edge.index) = std::make_pair(edge.vertices[0], edge.vertices[1]);
    }

    this->init(positions,
               radii,
               directions,
               densities,
               fixed,
               edges);
}

void Tree::iterate(double dt) {
    int m = this->rods.size();
    int n = this->constraints.size();

    // Clear lambda for each iteration
    std::vector<Vector<double, 6>> multipliers(n, Vector<double, 6>::Zero());

    std::vector<Vector3d> new_positions = this->iterate_predict_position(dt);
    std::vector<Quaterniond> new_orientations = this->iterate_predict_orientation(dt);

    // Record positions and orientations of fixed entities.
    std::vector<Vector3d> fixed_positions = new_positions;
    std::vector<Quaterniond> fixed_orientations = new_orientations;

    // this->generate_collision_constraints(new_positions);
    double error;
    double error_before = this->build_right_vector(new_positions,
                                                   new_orientations,
                                                   multipliers,
                                                   dt).cwiseAbs().maxCoeff();

    int iter = 0;
    do {
        // this->solve_gauss_seidel(new_positions, new_orientations, multipliers, dt);
        this->solve_direct_constraints(new_positions, new_orientations, multipliers, dt);
        error = this->build_right_vector(new_positions,
                                                       new_orientations,
                                                       multipliers,
                                                       dt).cwiseAbs().maxCoeff();
        // assert(error < error_before);
        ++iter;
    } while (error > 1e-10 && iter < 1);
    std::cout << "Error: " << error_before << " to " << error << std::endl;
    double error_after = this->build_right_vector(new_positions,
                                                   new_orientations,
                                                   multipliers,
                                                   dt).cwiseAbs().maxCoeff();

    // Update rods
    #pragma omp parallel for
    for (int j = 0; j < m; ++j) {
        if (rods[j].fixed) {
            rods[j].position = fixed_positions[j];
            rods[j].orientation = fixed_orientations[j];
        } else {
            // rods[j].velocity = (new_positions[j] - rods[j].position) / dt;
            rods[j].position = new_positions[j];

            // rods[j].angular_velocity = 2.0 * (rods[j].orientation.conjugate() * new_orientations[j]).vec() / dt;
            rods[j].orientation = new_orientations[j];
        }
    }
}

/**
 * @brief Apply initial prediction of position using previous velocity.
 * @param dt timestep
 * @return predicted positions
 */
std::vector<Vector3d> Tree::iterate_predict_position(double dt) const {
    int m = rods.size();
    std::vector<Vector3d> new_positions(m);
    for (const Rod &rod : this->rods) {
        Vector3d external_force = rod.mass() * this->gravity;
        new_positions[rod.index] = rod.position + dt * rod.velocity;
        // Apply force to non-fixed rods.
        if (!rod.fixed) {
            new_positions[rod.index] += dt * dt * (1.0 / rod.mass()) * external_force;
        }
        assert(new_positions[rod.index].allFinite());
    }
    return new_positions;
}

/**
 * @brief Apply initial prediction of orientation using half of ngular velocity.
 * @param dt timestep
 * @return predicted orientations
 */
std::vector<Quaterniond> Tree::iterate_predict_orientation(double dt) const {
    const int m = rods.size();
    std::vector<Quaterniond> new_orientations(m);
    #pragma omp parallel for
    for (int j = 0; j < m; ++j) {
        const Rod &rod = this->rods[j];
        const Quaterniond &q = rod.orientation;
        new_orientations[j] = q;
        new_orientations[j].coeffs() += dt * angular_velocity_to_quaternion(q) * rod.angular_velocity;
        new_orientations[j].normalize();
    }
    return new_orientations;
}



// Calculating the H matrix


void Tree::calculate_negative_jacobian(
        const std::vector<Quaterniond> &orientations,
        std::vector<Triplet<double>> &out) const {
    #pragma omp parallel for
    for (const Constraint &constraint : this->constraints) {
        int row = 6 * constraint.index;
        Darboux darboux(
                    orientations[constraint.rods[0]],
                orientations[constraint.rods[1]],
                constraint.median_length());

        // For each rod in the index, compute the Jacobian, and copy it to the
        // output.
        for (int constraint_rod = 0; constraint_rod < 2; ++constraint_rod) {
            // Calculate the jacobian within the appropriate block.
            int col = 6 * constraint.rods[constraint_rod];

            Vector3d p = (constraint_rod == 0) ? Vector3d(0, 0, 1) : Vector3d(0, 0, -1);
            p *= 0.5 * constraint.length[constraint_rod];

            Matrix<double, 6, 6> negative_jacobian = -constraint.jacobian(constraint_rod,
                                                                          p,
                                                                          orientations[constraint.rods[constraint_rod]],
                                                                          darboux);
            assert(negative_jacobian.allFinite());
            for (int r2 = 0; r2 < 6; ++r2) {
                for (int c2 = 0; c2 < 6; ++c2) {
                    if (r2 >= 3 && c2 < 3) {
                        assert(negative_jacobian(r2, c2) == 0);
                        continue;
                    }
                    #pragma omp critical
                    out.push_back(Triplet(row + r2, col + c2, negative_jacobian(r2, c2)));
                }
            }
        }
    }
}

void calculate_mass_matrix(const std::vector<Rod> &rods,
                           const std::vector<Quaterniond> &orientations,
                           std::vector<Triplet<double>> &out) {
    for (const Rod &rod : rods) {
        int start = 6 * rod.index;
        Matrix<double, 6, 6> rod_mass = rod.mass_matrix(orientations[rod.index]);
        assert(rod_mass.allFinite());
        // Copy mass matrix to the diagonal of out
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                if (row != col && (row < 3 || col < 3)) {
                    assert(rod_mass(row, col) == 0);
                    continue;
                }
                out.push_back(Triplet<double>(start + row, start + col, rod_mass(row, col)));
            }
        }
    }
}

void calculate_negative_compliance(
        const std::vector<Constraint> &constraints,
        double dt,
        std::vector<Triplet<double>> &out) {
    for (const Constraint &constraint : constraints) {
        int start = 6 * constraint.index;
        // Get the compliance matrix for the constraint
        Vector<double, 6> constraint_compliance = constraint.compliance_step(dt).diagonal();

        assert(constraint_compliance.allFinite());
        // Copy the compliance into the diagonal
        for (int j = 0; j < 6; ++j) {
            out.push_back(Eigen::Triplet<double>(start + j, start + j, -constraint_compliance[j]));
        }
    }
}

SparseMatrix<double> Tree::build_H_matrix(const std::vector<Vector3d> &new_positions,
                                          const std::vector<Quaterniond> &new_orientations,
                                          double dt) {
    const int rod_size = 6;
    const int rod_count = this->rods.size();
    const int constraint_size = 6;
    const int constraint_count = this->constraints.size();
    const int size = rod_size * rod_count + constraint_size * constraint_count;

    const int middle = rod_size * rod_count;

    std::vector<Triplet<double>> triplets;
    std::vector<Triplet<double>> jacobian_triplets;
    this->calculate_negative_jacobian(new_orientations, jacobian_triplets);
    // Adjust to the bottom left. Also copy to the top right.
    for (const Triplet<double> &triplet : jacobian_triplets) {
        assert(triplet.row() < constraint_size * constraint_count);
        assert(triplet.col() < rod_size * rod_count);
        triplets.push_back(Triplet<double>(triplet.row() + middle, triplet.col(), triplet.value()));
        triplets.push_back(Triplet<double>(triplet.col(), triplet.row() + middle, triplet.value()));
    }

    calculate_mass_matrix(this->rods, new_orientations, triplets);

    std::vector<Triplet<double>> negative_compliance_triplets;
    calculate_negative_compliance(this->constraints, dt, negative_compliance_triplets);
    for (const Triplet<double> &triplet : negative_compliance_triplets) {
        triplets.push_back(Triplet<double>(triplet.row() + middle, triplet.col() + middle, triplet.value()));
    }

    H.setFromTriplets(triplets.cbegin(), triplets.cend());
    return H;
}


VectorXd Tree::build_right_vector(const std::vector<Vector3d> &new_positions,
                                  const std::vector<Quaterniond> &new_orientations,
                                  const std::vector<Vector<double, 6>> &multipliers,
                                  double dt) const {
    const int rod_size = 6;
    const int rod_count = this->rods.size();
    const int constraint_size = 6;
    const int constraint_count = this->constraints.size();

    VectorXd negative_b(constraint_size * constraint_count);
    #pragma omp parallel for
    for (const Constraint &constraint : this->constraints) {
        negative_b.segment<constraint_size>(constraint_size * constraint.index) =
                constraint.evaluate(
                    new_positions,
                    new_orientations) +
                constraint.compliance_step(dt) * multipliers[constraint.index];
    }

    VectorXd right(rod_size * rod_count + constraint_size * constraint_count);
    right.head(rod_size * rod_count).setZero();
    right.tail(constraint_size * constraint_count) = negative_b;

    return right;
}

void Tree::solve_direct_constraints(std::vector<Vector3d> &new_positions,
                                    std::vector<Quaterniond> &new_orientations,
                                    std::vector<Vector<double, 6>> multipliers,
                                    double dt) {
    const int position_size = 3;
    const int orientation_size = 3;
    const int rod_size = position_size + orientation_size;
    const int rod_count = this->rods.size();
    const int constraint_size = 6;
    const int constraint_count = this->constraints.size();

    const SparseMatrix<double> H = this->build_H_matrix(new_positions, new_orientations, dt);
    const VectorXd right = this->build_right_vector(new_positions, new_orientations, multipliers, dt);
    assert(right.rows() == H.rows());
    assert(right.rows() == rod_size * rod_count + constraint_size * constraint_count);
    assert(H.coeffs().allFinite());

    if (!this->solver_prepared) {
        solver.analyzePattern(H);
        this->solver_prepared = true;
    }
    solver.factorize(H);
    const VectorXd solution = solver.solve(right);

    // The solution contains X and lambda terms.

    // Update x terms
    const auto dx = solution.head(rod_size * rod_count);
    #pragma omp parallel for
    for (int j = 0; j < rod_count; ++j) {
        // Skip fixed rods
        if (this->rods[j].fixed) {
            continue;
        }

        // Get start of this rod's values
        const int start = rod_size * j;
        // Copy position and orientation data
        new_positions[j] += dx.segment<position_size>(start);
        assert(new_positions[j].allFinite());

        const Matrix<double, 4, 3> G = angular_velocity_to_quaternion(new_orientations[j]);
        new_orientations[j].coeffs() += G * dx.segment<orientation_size>(start + position_size);
        assert(new_orientations[j].coeffs().allFinite());
        new_orientations[j].normalize();
    }

    // Update lambda terms
    const auto dlambda = solution.tail(constraint_size * constraint_count);
    #pragma omp parallel for
    for (int i = 0; i < constraint_count; ++i) {
        // Get start of this constraint's values
        int start = constraint_size * i;
        // Copy lambda
        multipliers[i] += dlambda.segment<constraint_size>(start);
        assert(multipliers[i].allFinite());
    }
}





void Tree::solve_gauss_seidel(std::vector<Vector3d> &new_positions,
                                    std::vector<Quaterniond> &new_orientations,
                                    std::vector<Vector<double, 6>> multipliers,
                                    double dt) {
    const int position_size = 3;
    const int orientation_size = 3;
    const int rod_size = position_size + orientation_size;
    const int rod_count = this->rods.size();
    const int constraint_size = 6;
    const int constraint_count = this->constraints.size();


    for (const Constraint &constraint : this->constraints) {
        Vector<double, 6> numerators = -constraint.evaluate(new_positions, new_orientations) - constraint.compliance_step(dt) * multipliers[constraint.index];

        for (int i = 0; i < 6; ++i) {
            double numerator = numerators[i];
        }
    }
}

}
