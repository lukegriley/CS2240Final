#include "rods/rod.h"

#include <iostream>
#include <stack>

#include <plant/plant.h>
#include "rods/util.h"
#include <omp.h>

using namespace Eigen;

namespace rod {

// Use the given vector in quaternion math.
Quaterniond as_quaternion(const Vector3d &vec) {
    return Quaterniond(0, vec[0], vec[1], vec[2]);
}

Vector3d Rod::direction(const Tree &tree) const  {
    return tree.particles[particles[1]].position - tree.particles[particles[0]].position;
}

void Rod::set_material_parameters(double youngs_modulus, double torsion_modulus) {
    double radius4 = this->radius * this->radius * this->radius * this->radius;
    double area_inertia = 0.25 * std::numbers::pi * radius4;

    DiagonalMatrix<double, 3> material = {
        area_inertia * youngs_modulus,
        area_inertia * youngs_modulus,
        2 * area_inertia * torsion_modulus
    };
    this->compliance = material.inverse();
}

Vector3d calculate_darboux(const Quaterniond &q1, const Quaterniond &q2) {
    return (q1.conjugate() * q2).vec();
}

void Tree::init_particles(std::vector<Vector3d> positions, std::vector<double> mass) {
    assert(positions.size() == mass.size());
    int n = positions.size();
    particles.resize(n);
    for (int i = 0; i < n; ++i) {
        particles[i].index = i;
        particles[i].position = positions[i];
        particles[i].velocity = Vector3d::Zero();
        particles[i].mass = mass[i];
        assert(mass[i] > 0);
        particles[i].fixed = false;
    }

}

void Tree::init_orientations(std::vector<std::pair<int, int>> rods, std::vector<double> radii) {
    assert(rods.size() == radii.size());
    using std::numbers::pi;
    int m = rods.size();
    this->rods.resize(m);
    lambda_shear.resize(m);
    lambda_twist.resize(m);
    // Compute particles and radius for each rod
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        rod.index = j;
        rod.particles[0] = std::get<0>(rods[j]);
        rod.particles[1] = std::get<1>(rods[j]);
        rod.radius = radii[j];
        rod.fixed = false;
        rod.rest_length = rod.direction(*this).norm();
        assert(rod.rest_length > 0);
    }

    // Compute parents for each rod
    int n = particles.size();
    std::vector<int> parents(n, -1);
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        assert(parents[rod.particles[1]] == -1);
        parents[rod.particles[1]] = j;
    }
    for (Rod &rod : this->rods) {
        rod.parent = parents[rod.particles[0]];
    }

    // Fix the rod's orientation at the root
    std::vector<bool> visited(m, false);
    std::vector<Quaterniond> orientations(m);
    std::stack<int> to_visit;
    for (const Rod &rod : this->rods) {
        to_visit.push(rod.index);
    }
    while (to_visit.size() > 0) {
        int next = to_visit.top();
        if (visited[next]) {
            to_visit.pop();
            continue;
        }
        Rod &rod = this->rods[next];
        int parent_index = rod.parent;

        // Set root, if the rod has no parent
        if (parent_index == -1) {
            orientations[next].setFromTwoVectors(Vector3d{0, 0, 1}, rod.direction(*this));
            to_visit.pop();
            visited[next] = true;
            continue;
        }

        // Visit parent first
        if (!visited[parent_index]) {
            to_visit.push(parent_index);
            continue;
        }

        const Rod &parent = this->rods[parent_index];
        // Initialize quaternion
        orientations[next] = Quaterniond::FromTwoVectors(parent.direction(*this), rod.direction(*this)) * orientations[parent_index];
        to_visit.pop();
        visited[next] = true;
    }
    for (int j = 0; j < m; ++j) {
        this->rods[j].orientation = orientations[j];
        this->rods[j].orientation.normalize();
        this->rods[j].initial_orientation = this->rods[j].orientation;
    }

    // Initialize angular velocity
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        // TODO: initializer angular velocity
        rod.angular_velocity = Vector3d::Zero();
        // TODO: initialize inertia (mass?)
    }


    // Compute simplified scalar inertia
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        // Masses of particles
        double masses[2] = {
            0.5 * particles[rod.particles[0]].mass,
            0.5 * particles[rod.particles[1]].mass,
        };
        double length = rod.direction(*this).norm();

        // Assume we rotate around the center of gravity.
        rod.inertia_s = masses[0] * masses[1] / (masses[0] + masses[1]) * length * length;
        // masses[0] * mean * mean +
        //        masses[1] * (length - mean) * (length - mean);
    }

    // Initialize inertia
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        rod.inertia = Matrix3d::Zero();
        rod.inertia.diagonal() = Vector3d(rod.inertia_s, rod.inertia_s, rod.inertia_s);
        rod.inverse_inertia = rod.inertia.inverse();
    }


    // For each particle, compute the Darboux vector.
    for (Rod &rod : this->rods) {
        if (rod.parent != -1) {
            const Rod &parent = this->rods[rod.parent];
            rod.initial_darboux = calculate_darboux(parent.orientation, rod.orientation);
        }
    }
}

void Tree::init_from_plant(const plant::Plant &plant, double density) {
    using std::numbers::pi;

    std::vector<Vector3d> positions(1 + plant.vertices.size());
    std::vector<double> masses(1 + plant.vertices.size(), 0);
    std::vector<std::pair<int, int>> rods(plant.vertices.size());
    std::vector<double> radii(plant.vertices.size());

    // Initialize first position to be zero
    positions[0] = Vector3d::Zero();

    for (const plant::Vertex &vertex : plant.vertices) {
        double volume = vertex.direction(plant).norm() * pi * vertex.radius * vertex.radius;
        double mass = volume * density;
        // Set the particle for the vertex's tail.
        positions.at(vertex.index + 1) = vertex.tail_position;
        masses.at(vertex.parent + 1) += mass / 2;
        masses.at(vertex.index + 1) += mass / 2;
        // Record the vertex as a rod.
        rods.at(vertex.index) = std::make_pair(vertex.parent + 1, vertex.index + 1);
        radii.at(vertex.index) = vertex.radius;
    }

    this->init_particles(positions, masses);
    this->init_orientations(rods, radii);

    // Fix the root
    this->particles[0].fixed = true;
    // Fix specified rods
    for (const plant::Vertex &vertex : plant.vertices) {
        if (vertex.fixed) {
            this->particles.at(vertex.parent + 1).fixed = true;
            this->particles.at(vertex.index + 1).fixed = true;
            this->rods.at(vertex.index).fixed = true;
        }
    }
    // Fix the root cylinder
    for (const plant::Vertex &vertex : plant.vertices) {
        if (vertex.parent == -1) {
            this->particles.at(vertex.parent + 1).fixed = true;
            this->particles.at(vertex.index + 1).fixed = true;
            this->rods.at(vertex.index).fixed = true;
        }
    }

}

void Tree::iterate(double dt) {
    int n = this->particles.size();
    int m = this->rods.size();

    // clear lambda for each iteration
    for(int i = 0; i < m; i++){
        lambda_shear[i] = Eigen::Vector3d::Zero();
        lambda_twist[i] = Eigen::Vector3d::Zero();
    }

    std::vector<Vector3d> new_velocities;
    this->iterate_predict_velocity(dt, new_velocities);
    std::vector<Vector3d> new_positions = this->iterate_predict_position(dt, new_velocities);
    std::vector<Vector3d> new_angular_velocities = this->iterate_predict_angular_velocity(dt);
    std::vector<Quaterniond> new_orientations = this->iterate_predict_orientation(dt, new_angular_velocities);

    // Record positions and orientations of fixed entities.
    std::vector<Vector3d> fixed_positions = new_positions;
    std::vector<Quaterniond> fixed_orientations = new_orientations;

    // this->generate_collision_constraints(new_positions);
    int iterations = 100;
    for (int iter = 0; iter < iterations; ++iter) {
        this->project_constraints(new_positions, new_orientations, dt);
    }

    // Update particles
    for (int i = 0; i < n; ++i) {
        if (particles[i].fixed) {
            particles[i].position = fixed_positions[i];
        } else {
            particles[i].velocity = (new_positions[i] - particles[i].position) / dt;
            particles[i].position = new_positions[i];
        }
    }
    // Update rods
    for (int j = 0; j < m; ++j) {
        if (rods[j].fixed) {
            rods[j].orientation = fixed_orientations[j];
        } else {
            rods[j].angular_velocity = 2.0 * (rods[j].orientation.conjugate() * new_orientations[j]).vec() / dt;
            rods[j].orientation = new_orientations[j];
        }
    }
}


void Tree::iterate_predict_velocity(double dt,
                                    std::vector<Vector3d> &new_velocities) const {
    int n = particles.size();
    new_velocities.resize(n);
    for (int i = 0; i < n; ++i) {
        new_velocities[i] = particles[i].velocity;
        if (!particles[i].fixed) {
            Vector3d external_force = this->particles[i].mass * this->gravity;
            new_velocities[i] += dt * particles[i].weight() * external_force;
        }
    }
}

/**
 * @brief Apply initial prediction of position using previous velocity.
 * @param dt timestep
 * @return predicted positions
 */
std::vector<Vector3d> Tree::iterate_predict_position(double dt, const std::vector<Vector3d> &new_velocities) const {
    int n = particles.size();
    std::vector<Vector3d> new_positions(n);
    for (int i = 0; i < n; ++i) {
        new_positions[i] = particles[i].position + dt * new_velocities[i];
    }
    return new_positions;
}

/**
 * @brief Apply initial prediction of angular velocity using external torque
 * @param dt timestep
 * @return predicted angular velocities
 */
std::vector<Vector3d> Tree::iterate_predict_angular_velocity(double dt) const {
    int m = rods.size();
    std::vector<Vector3d> new_angular_velocities(m);
    for (int j = 0; j < m; ++j) {
        const auto &omega = rods[j].angular_velocity;
        new_angular_velocities[j] = omega;
        if (!rods[j].fixed) {
            // Apply torque
            const Vector3d torque { 0, 0, 0 };
            new_angular_velocities[j] += dt * rods[j].inverse_inertia * (torque - omega.cross(rods[j].inertia * omega));
        }
    }
    return new_angular_velocities;
}

/**
 * @brief Apply initial prediction of orientation using half of ngular velocity.
 * @param dt timestep
 * @return predicted orientations
 */
std::vector<Quaterniond> Tree::iterate_predict_orientation(double dt,
                                                           const std::vector<Vector3d> &new_angular_velocities) const {
    int m = rods.size();
    std::vector<Quaterniond> new_orientations(m);
    for (int j = 0; j < m; ++j) {
        const Rod &rod = this->rods[j];
        const Quaterniond &q = rod.orientation;
        new_orientations[j] = q;
        Quaterniond angular_velocity = as_quaternion(new_angular_velocities[j]);
        new_orientations[j].coeffs() += 0.5 * dt * (q * angular_velocity).coeffs();
        new_orientations[j].normalize();
    }
    return new_orientations;
}

void Tree::generate_collision_constraints() {

}

void Tree::project_constraints(std::vector<Vector3d> &new_positions,
                               std::vector<Quaterniond> &new_orientations,
                               double dt) {
    this->project_stretch_shear_constraints(new_positions, new_orientations);
    this->project_bend_twist_constraints(new_orientations, dt);
}

void Tree::project_stretch_shear_constraints(std::vector<Vector3d> &new_positions,
                                             std::vector<Quaterniond> &new_orientations) {
    std::vector<Vector3d> dp(new_positions.size(), Vector3d::Zero());
    std::vector<Quaterniond> dq(new_orientations.size(), Quaterniond(0, 0, 0, 0));

    // Eigen::MatrixXd lambda = Eigen::MatrixXd::Zero(3, m);
    // Apply stretch-shear constraints
    int i = 0;
    for (const Rod &rod : this->rods) {
        double w1 = this->particles[rod.particles[0]].weight();
        double w2 = this->particles[rod.particles[1]].weight();
        double wq = rod.weight();
        const Vector3d &p1 = new_positions[rod.particles[0]];
        const Vector3d &p2 = new_positions[rod.particles[1]];
        double l = rod.rest_length;
        const Quaterniond &q = new_orientations[rod.index];
        Vector3d d = q * Vector3d(0, 0, 1);

        // --------------------------------------------------------------------------------------
        Quaterniond e3(0, 0, 0, 1);
        Quaterniond quat_z = e3 * q.conjugate();
        // Eigen::MatrixXd gradient_q = Eigen::MatrixXd::Zero(3, 4);
        /*
         * for q = w + xi + yj + zk
         * S(v) =
         *          [0,  -z, y
         *           z,  0 , -x
         *           -y, x , 0]
         * so -s(v) =
         *          [0,  z, -y
         *           -z, 0 , x
         *           y,  x , 0]
         */
        // gradient_q(0,0) = quaternion_z.x();
        // gradient_q(1,0) = quaternion_z.y();
        // gradient_q(2,0) = quaternion_z.z();
        // gradient_q(0,1) = quaternion_z.w();
        // gradient_q(1,2) = quaternion_z.w();
        // gradient_q(2,3) = quaternion_z.w();
        // gradient_q(0,2) = quaternion_z.z();
        // gradient_q(0,3) = -quaternion_z.y();
        // gradient_q(1,1) = -quaternion_z.z();
        // gradient_q(1,3) = quaternion_z.x();
        // gradient_q(2,1) = quaternion_z.y();
        // gradient_q(2,2) = -quaternion_z.x();
        // gradient_q *= 2;

        // MatrixXd prod1 = gradient_q * gradient_q.transpose();
        // Vector3d q_diag(prod1(0,0), prod1(1,1), prod1(2,2));
        // // double denominator = w1 + w2 + 4 * wq * l * l;  //denominator shoud be vec3, wq(2) is different from wq(0) and wq(1)
        // Vector4d prod2 = gradient_q.transpose() * diff;
        // Quaterniond g_cst_st(prod2(0), prod2(1), prod2(2), prod2(3));
        // Quaterniond g_cst_st(w, x, y, z);
        // dq[rod.index].coeffs() += (wq * l * l) / denominator * g_cst_st.coeffs();

        // ------------------------------------------------------------------------------------------------------
        //alpha_tilde = alpha_inverse / (t^2); the paper suggests alpha_inverse [1e-8, 1e-10], our delta_t = 1e-3
        const double alpha_tilde = 0.01;

        // Vector3d intertias(wq, wq, 2*wq);
        Vector3d shear_constraint_eval = 1 / l * (p2 - p1) - d;

        double common_in_denom = w1 + w2 + alpha_tilde * l * l + 4 * wq * l * l;
        // Vector3d denominators(common_in_denom, common_in_denom, common_in_denom + 4 * wq * l * l);  //J = I1 + I2 = 2I1 in our case
        // Vector3d denominators(common_in_denom, common_in_denom, common_in_denom);  //J = I1 + I2 = 2I1 in our case
        Vector3d numerators = l * l * shear_constraint_eval + l * l * alpha_tilde * lambda_shear[i];
        // Vector3d delta_lambda = numerators.array() / denominators.array();
        Vector3d delta_lambda = numerators / common_in_denom;
        // update lambda
        lambda_shear[i] += delta_lambda;
        // update orientation
        // dc1/dq0 * d_lambda + dc2/dq0 * d_lambda + dc3/dq0 * d_lambda
        double w = quat_z.x() * delta_lambda.x() + quat_z.y() * delta_lambda.y() + quat_z.z() * delta_lambda.z();
        double x = quat_z.w() * delta_lambda.x() - quat_z.z() * delta_lambda.y() + quat_z.y() * delta_lambda.z();
        double y = quat_z.z() * delta_lambda.x() + quat_z.w() * delta_lambda.y() - quat_z.x() * delta_lambda.z();
        double z = -quat_z.y() * delta_lambda.x() + quat_z.x() * delta_lambda.y() + quat_z.w() * delta_lambda.z();
        // Quaterniond quat_dq(wq * w, wq * x, wq * y, wq * z);
        // dq[rod.index].coeffs() += quat_dq.coeffs();
        Quaterniond quat_dq(w, x, y, z);
        dq[rod.index].coeffs() += wq * 2 * quat_dq.coeffs();
        // // update position
        dp[rod.particles[0]] += w1 / l * delta_lambda;
        dp[rod.particles[1]] -= w2 / l * delta_lambda;
        assert(dp[rod.particles[0]].allFinite());
        assert(dp[rod.particles[1]].allFinite());
        i++;
    }

    for (int i = 0; i < this->particles.size(); ++i) {
        new_positions[i] += dp[i];
    }
    for (const Rod &rod : rods) {
        assert(dq[rod.index].coeffs().allFinite());
        new_orientations[rod.index].coeffs() += dq[rod.index].coeffs();
        new_orientations[rod.index].normalize();
    }
}


void Tree::project_bend_twist_constraints(std::vector<Quaterniond> &new_orientations,
                                          double dt) {
    // Apply bend-twist constraints
    // Do this over multiple iterations for stability
    // We apply these in an interleaving order.
    for (int j = 0; j < this->rods.size(); ++j) {
        // Note: the +1 and -1 are there to skip over the first fixed rod.
        // In the future, we would want to generalize this to trees.
        int rod_index = interleave(this->rods.size(), j);
        // rod_index = std::rand() % this->rods.size();
        const Rod &rod = this->rods.at(rod_index);
        if (rod.parent == -1) {
            continue;
        }
        Quaterniond q = new_orientations[rod.parent];
        Quaterniond u = new_orientations[rod.index];

        int num_steps = num_bend_twist_steps;
        for (int i = 0; i < num_steps; ++i) {
            Quaterniond dq, du;
            std::tie(dq, du) = this->project_bend_twist_constraint(rod, q, u, dt, rod.index);

            q = (q.coeffs() + 1.0 / num_steps * dq.coeffs()).normalized();
            u = (u.coeffs() + 1.0 / num_steps * du.coeffs()).normalized();

        }

        new_orientations[rod.parent] = q;
        new_orientations[rod.index] = u;
    }

}


std::pair<Quaterniond, Quaterniond> Tree::project_bend_twist_constraint(
        const Rod &rod,
        const Quaterniond &q,
        const Quaterniond &u,
        double dt,
        int iter) {

    assert(rod.parent != -1);
    const Rod &parent = this->rods.at(rod.parent);

    Vector3d darboux = calculate_darboux(q, u);
    Vector3d initial_darboux = rod.initial_darboux;
    // We suppose q is the parent, u is this rod
    double s = (darboux.dot(initial_darboux) > 0) ? 1 : -1;
    double wq = parent.weight();
    double wu = rod.weight();

    Eigen::Quaterniond dq, du;

    Quaterniond darboux_diff = as_quaternion(darboux - s * initial_darboux);

    // --------------------------------------------------------
    // Eigen::MatrixXd gradient_q = Eigen::MatrixXd::Zero(3, 4);
    // gradient_q(0,0) = -u.x();
    // gradient_q(1,0) = -u.y();
    // gradient_q(2,0) = -u.z();
    // gradient_q(0,1) = u.w();
    // gradient_q(1,2) = u.w();
    // gradient_q(2,3) = u.w();
    // gradient_q(0,2) = u.z();
    // gradient_q(0,3) = -u.y();
    // gradient_q(1,1) = -u.z();
    // gradient_q(1,3) = u.x();
    // gradient_q(2,1) = u.y();
    // gradient_q(2,2) = -u.x();

    // Eigen::MatrixXd gradient_u = Eigen::MatrixXd::Zero(3, 4);
    // gradient_u(0,0) = -q.x();
    // gradient_u(1,0) = -q.y();
    // gradient_u(2,0) = -q.z();
    // gradient_u(0,1) = q.w();
    // gradient_u(1,2) = q.w();
    // gradient_u(2,3) = q.w();
    // gradient_u(0,2) = q.z();
    // gradient_u(0,3) = -q.y();
    // gradient_u(1,1) = -q.z();
    // gradient_u(1,3) = q.x();
    // gradient_u(2,1) = q.y();
    // gradient_u(2,2) = -q.x();
    // gradient_u *= -1;

    // Vector4d prod_q = gradient_q.transpose() * (darboux - s * initial_darboux);
    // Quaterniond g_cst_st_q(prod_q(0),prod_q(1),prod_q(2),prod_q(3));

    // Vector4d prod_u = gradient_u.transpose() * (darboux - s * initial_darboux);
    // Quaterniond g_cst_st_u(prod_u(0),prod_u(1),prod_u(2),prod_u(3));

    // dq.coeffs() = wq / (wq + wu) * (g_cst_st_q).coeffs();
    // du.coeffs() = wu / (wq + wu) * (g_cst_st_u).coeffs();

    // --------------------------------------------------------
    Vector3d twist_constraint_eval = darboux - s * initial_darboux;
    Vector3d c = twist_constraint_eval;

    // w := EI1 + EI2 + GJ
    // denominators_pbd :=
    //          (pho * l_q * I1_q + pho * l_u * I1_u,
    //           pho * l_q * I2_q + pho * l_u * I2_u,
    //           pho * l_q * J_q  + pho * l_u * J_u)
    // alpha := (EI1, EI2, GJ)
    // alpha_tilde = alpha / (dt * dt)

    Vector3d alpha_tilde = rod.compliance.diagonal() / (dt * dt);
    Vector3d denominator_pbd(wq + wu, wq + wu, wq + wu);
    Vector3d denominator = denominator_pbd + alpha_tilde;
    Vector3d numerator = twist_constraint_eval - alpha_tilde.cwiseProduct(lambda_twist[rod.index]);
    Vector3d delta_lambda = numerator.cwiseQuotient(denominator);
    // update lambda
    lambda_twist[rod.index] += delta_lambda;
    // Similar to shear, dc1/dq0 * d_lambda + dc2/dq0 * d_lambda + dc3/dq0 * d_lambda for q
    double q_w = -u.x() * delta_lambda.x() -u.y() * delta_lambda.y() - u.z() * delta_lambda.z();
    double q_x = u.w() * delta_lambda.x() - u.z() * delta_lambda.y() + u.y() * delta_lambda.z();
    double q_y = u.z() * delta_lambda.x() + u.w() * delta_lambda.y() - u.x() * delta_lambda.z();
    double q_z = -u.y() * delta_lambda.x() + u.x() * delta_lambda.y() + u.w() * delta_lambda.z();
    Quaterniond quat_dq(q_w, q_x, q_y, q_z);
    dq.coeffs() = wq * quat_dq.coeffs();
    // Similar to dq, dc1/dq0 * d_lambda + dc2/dq0 * d_lambda + dc3/dq0 * d_lambda for q
    double u_w = -q.x() * delta_lambda.x() -q.y() * delta_lambda.y() - q.z() * delta_lambda.z();
    double u_x = q.w() * delta_lambda.x() - q.z() * delta_lambda.y() + q.y() * delta_lambda.z();
    double u_y = q.z() * delta_lambda.x() + q.w() * delta_lambda.y() - q.x() * delta_lambda.z();
    double u_z = -q.y() * delta_lambda.x() + q.x() * delta_lambda.y() + q.w() * delta_lambda.z();
    Quaterniond quat_du(u_w, u_x, u_y, u_z);
    du.coeffs() = -wu * quat_du.coeffs();

    assert(dq.coeffs().allFinite());
    assert(dq.coeffs().allFinite());
    return std::make_pair(dq, du);
}

}
