#pragma once

#include <Eigen/Dense>
#include <plant/plant.h>

namespace rod {

struct Particle;
struct Rod;
struct Tree;

struct Particle {
    int index;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d initial_darboux;
    double mass;
    inline double weight() const { return 1.0 / mass; }
    bool fixed;
};

struct Rod {
    int index = -1;
    double radius;
    double rest_length;
    Eigen::Quaterniond initial_orientation;
    Eigen::Quaterniond orientation;
    Eigen::Vector3d angular_velocity;
    Eigen::Matrix3d inertia;
    Eigen::Matrix3d inverse_inertia;
    double inertia_s;
    inline double weight() const { return 1.0 / inertia_s; }

    Eigen::Vector3d initial_darboux;

    bool fixed = false;
    int parent = -1;
    int particles[2];

    Eigen::Vector3d direction(const Tree &tree) const;
};

struct Tree {
    std::vector<Particle> particles;
    std::vector<Rod> rods;
    std::vector<Eigen::Vector3d> lambda_shear;
    std::vector<Eigen::Vector3d> lambda_twist;
    void init_particles(std::vector<Eigen::Vector3d> positions, std::vector<double> mass);
    void init_orientations(std::vector<std::pair<int, int>> rods, std::vector<double> radii);
    void init_from_plant(const plant::Plant &plant, double density);
    void iterate(double dt);

    int num_bend_twist_steps = 20;
    Eigen::Vector3d gravity {0, 0, 0};

private:
    void iterate_predict_velocity(double dt,
                                  std::vector<Eigen::Vector3d> &new_velocities) const;

    /**
     * @brief Apply initial prediction of position using previous velocity.
     * @param dt timestep
     * @param new_velocities predicted velocities
     * @return predicted positions
     */
    std::vector<Eigen::Vector3d> iterate_predict_position(double dt, const std::vector<Eigen::Vector3d> &new_velocities) const;

    /**
     * @brief Apply initial prediction of angular velocity using external torque
     * @param dt timestep
     * @return predicted angular velocities
     */
    std::vector<Eigen::Vector3d> iterate_predict_angular_velocity(double dt) const;

    /**
     * @brief Apply initial prediction of orientation using half of ngular velocity.
     * @param dt timestep
     * @return predicted orientations
     */
    std::vector<Eigen::Quaterniond> iterate_predict_orientation(double dt,
                                                                const std::vector<Eigen::Vector3d> &new_angular_velocities) const;

    void generate_collision_constraints();
    void project_constraints(std::vector<Eigen::Vector3d> &new_positions,
                             std::vector<Eigen::Quaterniond> &new_orientations);
    void project_stretch_shear_constraints(std::vector<Eigen::Vector3d> &new_positions,
                             std::vector<Eigen::Quaterniond> &new_orientations);
    void project_bend_twist_constraints(std::vector<Eigen::Quaterniond> &new_orientations);
    std::pair<Eigen::Quaterniond, Eigen::Quaterniond> project_bend_twist_constraint(
            const Rod &rod,
            const Eigen::Quaterniond &q, const Eigen::Quaterniond &u,
            int iter);
};

}
