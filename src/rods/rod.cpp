#include "rod.h"

#include <iostream>

using namespace Eigen;

// Use the given vector in quaternion math.
Quaterniond as_quaternion(const Vector3d &vec) {
    return Quaterniond(0, vec[0], vec[1], vec[2]);
}

Vector3d Rod::direction(const Tree &tree) const  {
    return tree.particles[particles[1]].position - tree.particles[particles[0]].position;
}

Vector3d Rod::darboux(const Tree &tree, const std::vector<Quaterniond> &orientations) const {
    if (this->parent == -1) {
        throw std::runtime_error("no parent");
    }

    const Rod &parent = tree.rods[this->parent];
    double l = 0.5 * (this->initial_direction.norm() + parent.initial_direction.norm());
    const Quaterniond &parent_orientation = orientations[this->parent];
    const Quaterniond &this_orientation = orientations[this->index];

    return 2 / l * (parent_orientation.conjugate() * this_orientation).vec();
}
Vector3d Rod::initial_darboux(const Tree &tree) const {
    if (this->parent == -1) {
        throw std::runtime_error("no parent");
    }
    const Rod &parent = tree.rods[this->parent];
    double l = 0.5 * (this->initial_direction.norm() + parent.initial_direction.norm());
    const Quaterniond &parent_orientation = parent.initial_orientation;
    const Quaterniond &this_orientation = this->initial_orientation;

    return 2 / l * (parent_orientation.conjugate() * this_orientation).vec();
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
        particles[i].fixed = false;
    }
}

void Tree::init_orientations(std::vector<std::pair<int, int>> rods, std::vector<double> radii) {
    assert(rods.size() == radii.size());
    using std::numbers::pi;
    int m = rods.size();
    this->rods.resize(m);
    // Compute particles and radius for each rod
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        rod.index = j;
        rod.particles[0] = std::get<0>(rods[j]);
        rod.particles[1] = std::get<1>(rods[j]);
        rod.radius = radii[j];
        rod.fixed = false;
    }


    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        rod.initial_direction = rod.direction(*this);
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
    this->rods[0].orientation.setFromTwoVectors(Vector3d{0, 0, 1}, this->rods[0].direction(*this));
    for (int j = 1; j < m; ++j) {
        // Initialize quaternions from position
        // TODO: ensure that the parent is initialized.
        Rod &rod = this->rods[j];
        const Rod &parent = this->rods[rod.parent];
        // Find the transformation nearest to the parent.
        rod.orientation = Quaterniond::FromTwoVectors(parent.direction(*this), rod.direction(*this)) * parent.orientation;
        rod.orientation.normalize();
        rod.initial_orientation = rod.orientation;
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
        double mass = masses[0] + masses[1];
        double length = rod.direction(*this).norm();
        double mean = masses[1] / mass * length;

        // Assume we rotate around the center of gravity.
        rod.inertia_s = masses[0] * masses[1] / (masses[0] + masses[1]) * length * length;
        // masses[0] * mean * mean +
        //        masses[1] * (length - mean) * (length - mean);
    }

    // Initialize inertia
    for (int j = 0; j < m; ++j) {
        Rod &rod = this->rods[j];
        double mass = 0.5 * (particles[rod.particles[0]].mass + particles[rod.particles[1]].mass);
        double radius = 1;
        double radius4 = radius * radius * radius * radius;
        rod.inertia = Matrix3d::Zero();
        rod.inertia.diagonal() = Vector3d(rod.inertia_s, rod.inertia_s, rod.inertia_s);
        rod.inverse_inertia = rod.inertia.inverse();
    }


    // For each particle, compute the Darboux vector.
}

void Tree::iterate(double dt) {
    int n = this->particles.size();
    int m = this->rods.size();

    std::vector<Vector3d> new_velocities;
    this->iterate_predict_velocity(dt, new_velocities);
    std::vector<Vector3d> new_positions = this->iterate_predict_position(dt, new_velocities);
    std::vector<Vector3d> new_angular_velocities = this->iterate_predict_angular_velocity(dt);
    std::vector<Quaterniond> new_orientations = this->iterate_predict_orientation(dt, new_angular_velocities);

    // this->generate_collision_constraints(new_positions);
    int iterations = 1000;
    for (int iter = 0; iter < iterations; ++iter) {
        this->project_constraints(new_positions, new_orientations);
    }

    // Update particles
    for (int i = 0; i < n; ++i) {
        if (particles[i].fixed)
            continue;
        particles[i].velocity = (new_positions[i] - particles[i].position) / dt;
        particles[i].position = new_positions[i];
    }
    // Update rods
    for (int j = 0; j < m; ++j) {
        if (rods[j].fixed)
            continue;
        rods[j].angular_velocity = 2.0 * (rods[j].orientation.conjugate() * new_orientations[j]).vec() / dt;
        rods[j].orientation = new_orientations[j];
    }
}


void Tree::iterate_predict_velocity(double dt,
                                    std::vector<Vector3d> &new_velocities) const {
    int n = particles.size();
    new_velocities.resize(n);
    for (int i = 0; i < n; ++i) {
        new_velocities[i] = particles[i].velocity;
        if (!particles[i].fixed) {
            Vector3d external_force = this->particles[i].mass * Vector3d(0.0, 0.0, -1);
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
        new_positions[i] = particles[i].position;
        if (!particles[i].fixed) {
            new_positions[i] += dt * new_velocities[i];
        }
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
        if (rods[j].fixed) {
            new_angular_velocities[j].setZero();
        } else {
            const auto &omega = rods[j].angular_velocity;
            const Vector3d torque { 0, 0, 0 };
            new_angular_velocities[j] = omega +
                    dt * rods[j].inverse_inertia * (torque - omega.cross(rods[j].inertia * omega));
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
        if (!rod.fixed) {
            Quaterniond angular_velocity = as_quaternion(new_angular_velocities[j]);
            new_orientations[j].coeffs() += 0.5 * dt * (q * angular_velocity).coeffs();
            new_orientations[j] = new_orientations[j].normalized();
        }
    }
    return new_orientations;
}

void Tree::generate_collision_constraints() {

}

void Tree::project_constraints(std::vector<Vector3d> &new_positions,
                               std::vector<Quaterniond> &new_orientations) const {
    this->project_stretch_shear_constraints(new_positions, new_orientations);
    this->project_bend_twist_constraints(new_positions, new_orientations);
}

void Tree::project_stretch_shear_constraints(std::vector<Vector3d> &new_positions,
                                             std::vector<Quaterniond> &new_orientations) const {
    std::vector<Vector3d> dp(new_positions.size(), Vector3d::Zero());
    std::vector<Quaterniond> dq(new_orientations.size(), Quaterniond(0, 0, 0, 0));

    // Apply stretch-shear constraints
    for (const Rod &rod : this->rods) {
        double w1 = this->particles[rod.particles[0]].weight();
        double w2 = this->particles[rod.particles[1]].weight();
        double wq = rod.weight();
        const Vector3d &p1 = new_positions[rod.particles[0]];
        const Vector3d &p2 = new_positions[rod.particles[1]];
        double l = rod.initial_direction.norm();
        const Quaterniond &q = new_orientations[rod.index];
        // TODO: check if this should be initial_orientation, old positions, or new positions.
        // Vector3d d = (p2 - p1).normalized(); // Looks liek a pendulum
        Vector3d d = q * Vector3d(0, 0, 1);
        Quaterniond conj_e3(0, 0, 0, -1);

        double denominator = w1 + w2 + 4 * wq * l * l;
        Vector3d diff = 1 / l * (p2 - p1) - d;
        Quaterniond diffq = as_quaternion(diff);
        dp[rod.particles[0]] += (w1 * l) / denominator * diff;
        dp[rod.particles[1]] -= (w2 * l) / denominator * diff;
        dq[rod.index].coeffs() += (wq * l * l) / denominator * (diffq * q * conj_e3).coeffs();
    }

    for (int i = 0; i < this->particles.size(); ++i) {
        new_positions[i] += dp[i];
    }
    for (const Rod &rod : rods) {
        new_orientations[rod.index].coeffs() += dq[rod.index].coeffs();
        new_orientations[rod.index].normalize();
    }
}


void Tree::project_bend_twist_constraints(std::vector<Vector3d> &new_positions,
                                             std::vector<Quaterniond> &new_orientations) const {
    std::vector<Vector3d> dp(new_positions.size(), Vector3d::Zero());
    std::vector<Quaterniond> dq(new_orientations.size(), Quaterniond(0, 0, 0, 0));

    // Apply bend-twist constraints
    for (const Rod &rod : this->rods) {
        if (rod.parent == -1) {
            continue;
        }
        const Rod &parent = this->rods.at(rod.parent);

        Vector3d darboux = rod.darboux(*this, new_orientations);
        Vector3d initial_darboux = rod.initial_darboux(*this);
        // We suppose q is the parent, u is this rod
        const Quaterniond &q = new_orientations.at(rod.parent);
        const Quaterniond &u = new_orientations.at(rod.index);
        double s = (darboux.dot(initial_darboux) > 0) ? 1 : -1;
        double wq = parent.weight();
        double wu = rod.weight();
        Quaterniond darboux_diff = as_quaternion(darboux - s * initial_darboux);
        // Hack fix
        darboux_diff.coeffs() = 0.01 * darboux_diff.coeffs();
        dq.at(rod.parent).coeffs() += wq / (wq + wu) * (u * darboux_diff).coeffs();
        dq.at(rod.index).coeffs() -= wu / (wq + wu) * (q * darboux_diff).coeffs();
    }

    for (int i = 0; i < this->particles.size(); ++i) {
        new_positions[i] += dp[i];
    }
    for (const Rod &rod : rods) {
        new_orientations[rod.index].coeffs() += dq[rod.index].coeffs();
        new_orientations[rod.index].normalize();
    }
}
