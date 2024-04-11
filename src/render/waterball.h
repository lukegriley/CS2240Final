#include "graphics/shader.h"
#include "graphics/shape.h"
#include "plant/plant.h"

struct WaterBallPlantRenderer {
    WaterBallPlantRenderer() = default;
    WaterBallPlantRenderer(const Plant &plant);
    WaterBallPlantRenderer(const WaterBallPlantRenderer &renderer) = delete;

    void init(const Plant &plant);
    void update(const Plant &plant);
    void render(Shader *shader);
private:
    int vertex_count;
    std::vector<Eigen::Affine3d> transformations;
    std::vector<Eigen::Vector3f> colors;
    Shape m_ball;
    float uniform_scale = 10.f;
};
