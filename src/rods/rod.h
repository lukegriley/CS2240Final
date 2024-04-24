#pragma once

#include <Eigen/Dense>

struct Particle;
struct Rod;
struct Tree;

struct Constraint {
    virtual double operator()(
            const Tree &tree,
            const std::vector<Eigen::Vector3d> &new_positions,
            const std::vector<Eigen::Vector3d> &new_orientations) const = 0;
    virtual std::vector<std::pair<int, Eigen::Vector3d>>
    position_gradient(
            int index,
            const Tree &tree,
            const std::vector<Eigen::Vector3d> &new_positions,
            const std::vector<Eigen::Vector3d> &new_orientations) const = 0;
    virtual std::vector<std::pair<int, Eigen::Quaterniond>>
    orientation_gradient(
            int index,
            const Tree &tree,
            const std::vector<Eigen::Vector3d> &new_positions,
            const std::vector<Eigen::Vector3d> &new_orientations) const = 0;
};


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
    int index;
    double radius;
    Eigen::Vector3d initial_direction;
    Eigen::Quaterniond orientation;
    Eigen::Vector3d angular_velocity;
    Eigen::Matrix3d inertia;
    Eigen::Matrix3d inverse_inertia;
    double inertia_s;
    inline double weight() const { return 1.0 / inertia_s; }

    bool fixed;
    int parent;
    int particles[2];

    Eigen::Vector3d direction(const Tree &tree) const;

};

struct Tree {
    std::vector<Particle> particles;
    std::vector<Rod> rods;

    void init_particles(std::vector<Eigen::Vector3d> positions, std::vector<double> mass);
    void init_orientations(std::vector<std::pair<int, int>> rods, std::vector<double> radii);
    void iterate(double dt);

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
};

