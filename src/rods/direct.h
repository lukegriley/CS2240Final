#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <plant/plant.h>

namespace rod::direct {

struct Rod;
struct Tree;

struct Rod {
    int index;
    double radius;
    double length;
    double density;
    bool fixed;

    int parent = -1;
    Eigen::Vector3d position;
    Eigen::Quaterniond orientation;
    Eigen::Vector3d velocity;
    Eigen::Vector3d angular_velocity;

    Rod() = default;
    Rod(int index,
        const Eigen::Vector3d &position,
        double radius,
        double length,
        double density,
        bool fixed);

    double volume() const;
    double mass() const;
    Eigen::Vector3d direction() const;
    double area_inertia() const;
    Eigen::Matrix<double, 6, 6> mass_matrix(const Eigen::Quaterniond &orientation) const;

    void set_material_parameters(double youngs_modulus, double torsion_modulus);
};

struct Darboux {
    std::array<Eigen::Quaterniond, 2> q;
    double length;

    Darboux(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2, double length);
    Darboux(const Darboux &other) = delete;
    Eigen::Vector3d vec() const;
    Eigen::Matrix<double, 3, 4> differentiate(int rod_index) const;
};


struct Constraint {
    int index;
    std::array<int, 2> rods;
    std::array<double, 2> length;
    std::array<double, 2> radii;

    Eigen::Vector3d initial_darboux;
    Eigen::DiagonalMatrix<double, 6> compliance; // alpha

    Constraint() = default;
    Constraint(int index, const Rod &rod1, const Rod &rod2);
    double median_length() const;
    void set_material_parameters(double youngs_modulus, double torsion_modulus);
    Eigen::DiagonalMatrix<double, 6> compliance_step(double dt) const; // tilde alpha

    Eigen::Vector<double, 6> evaluate(
            const std::vector<Eigen::Vector3d> &positions,
            const std::vector<Eigen::Quaterniond> &orientations) const;
    Eigen::Matrix<double, 6, 6> jacobian(
            int rod_index,
            const Eigen::Vector3d &p,
            const Eigen::Quaterniond &q,
            const Darboux &darboux) const;
};


struct Tree {
    std::vector<Rod> rods;
    std::vector<Constraint> constraints;

    void init(const std::vector<Eigen::Vector3d> &positions,
              const std::vector<double> &radii,
              const std::vector<Eigen::Vector3d> &directions,
              const std::vector<double> &masses,
              const std::vector<bool> &fixed,
              const std::vector<std::pair<int, int>> &edges);
    void init_from_plant(const plant::Plant &plant, double density);
    void compute_rest_orientation(const std::vector<Eigen::Vector3d> &directions,
                                  const std::vector<int> &parents);

    void iterate(double dt);

    Eigen::Vector3d gravity {0, 0, 0};

private:
    /**
     * @brief Apply initial prediction of position using previous velocity.
     * @param dt timestep
     * @param new_velocities predicted velocities
     * @return predicted positions
     */
    std::vector<Eigen::Vector3d> iterate_predict_position(double dt) const;

    /**
     * @brief Apply initial prediction of orientation using half of ngular velocity.
     * @param dt timestep
     * @return predicted orientations
     */
    std::vector<Eigen::Quaterniond> iterate_predict_orientation(double dt) const;

    void generate_collision_constraints();

    void solve_gauss_seidel(
            std::vector<Eigen::Vector3d> &new_positions,
            std::vector<Eigen::Quaterniond> &new_orientations,
            std::vector<Eigen::Vector<double, 6>> multipliers,
            double dt);
    void solve_direct_constraints(
            std::vector<Eigen::Vector3d> &new_positions,
            std::vector<Eigen::Quaterniond> &new_orientations,
            std::vector<Eigen::Vector<double, 6>> multipliers,
            double dt);

    Eigen::SparseMatrix<double> H;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool solver_prepared = false;
    Eigen::SparseMatrix<double> build_H_matrix(
            const std::vector<Eigen::Vector3d> &new_positions,
            const std::vector<Eigen::Quaterniond> &new_orientations,
            double dt);
    Eigen::VectorXd build_right_vector(
            const std::vector<Eigen::Vector3d> &new_positions,
            const std::vector<Eigen::Quaterniond> &new_orientations,
            const std::vector<Eigen::Vector<double, 6>> &multipliers,
            double dt) const;

    void calculate_negative_jacobian(
            const std::vector<Eigen::Quaterniond> &orientations,
            std::vector<Eigen::Triplet<double>> &out) const;
};

}
