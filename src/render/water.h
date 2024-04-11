#include "plant/plant.h"
#include "render/line.h"
#include "render/waterball.h"

struct WaterPlantRenderer {
    WaterPlantRenderer() = default;
    WaterPlantRenderer(const Plant &plant);
    WaterPlantRenderer(const WaterBallPlantRenderer &renderer) = delete;

    void init(const Plant &plant);
    void update(const Plant &plant);
    void render(Shader *shader);
private:
    LinePlantRenderer stem_renderer;
    WaterBallPlantRenderer vertex_renderer;
};
