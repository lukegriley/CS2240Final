#include "config.h"

// Define conversion rules

template<class T>
struct QVariantType {
    constexpr static std::string name = "unknown type";
    static inline T convert(const QVariant &) {
        throw std::runtime_error("Failed to convert type " + QVariantType<T>::name);
    }
};

template<>
struct QVariantType<std::string> {
    constexpr static std::string name = "string";
    static inline std::string convert(const QVariant &variant) {
        return variant.toString().toStdString();
    }
};

template<>
struct QVariantType<int> {
    constexpr static std::string name = "int";
    static inline int convert(const QVariant &variant) {
        return variant.toInt();
    }
};

template<>
struct QVariantType<double> {
    constexpr static std::string name = "double";
    static inline double convert(const QVariant &variant) {
        return variant.toDouble();
    }
};


/**
 * @brief check_key Check that the given key is defined
 */
template<class T>
void check_key(const QSettings &settings, const std::string &key) {
    if (!settings.contains(key)) {
        throw std::runtime_error("missing key " + key +
                                 " of type " + QVariantType<T>::name);
    }
}

// Convert a key to a certain value
template<class T>
T get_setting(const QSettings &settings, const std::string &key) {
    check_key<T>(settings, key);
    QVariant variantValue = settings.value(key);
    return QVariantType<T>::convert(variantValue);
}


void Config::init(const QSettings &settings) {
    this->io.plant_file = get_setting<std::string>(settings, "IO/plant");
    this->io.output_directory = get_setting<std::string>(settings, "IO/output");

    this->simulation.frames = get_setting<int>(settings, "Simulation/frames");

    this->water.dynamic_viscosity = get_setting<double>(settings, "Water/dynamicViscosity");
    this->water.uniform_loss_rate = get_setting<double>(settings, "Water/uniformLossRate");
    this->water.time_between_frames = get_setting<double>(settings, "Water/timeBetweenFrames");
    this->water.integration_time_step = get_setting<double>(settings, "Water/integrationTimeStep");
    this->water.elastic_steepness = get_setting<double>(settings, "Water/elasticSteepness");
    this->water.elastic_threshold = get_setting<double>(settings, "Water/elasticThreshold");

    this->rods.density = get_setting<double>(settings, "Rods/density");
    this->rods.gravity = get_setting<double>(settings, "Rods/gravity");
    this->rods.iterations_per_frame = get_setting<int>(settings, "Rods/iterationsPerFrame");
    this->rods.iteration_time_step = get_setting<double>(settings, "Rods/iterationTimeStep");
    this->rods.num_constraint_iterations = get_setting<int>(settings, "Rods/numConstraintIterations");

    this->rods.min_youngs_modulus = get_setting<double>(settings, "Rods/minYoungsModulus");
    this->rods.max_youngs_modulus = get_setting<double>(settings, "Rods/maxYoungsModulus");
    this->rods.torsion_modulus = get_setting<double>(settings, "Rods/torsionModulus");
}
