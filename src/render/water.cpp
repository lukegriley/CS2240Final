#include "water.h"

WaterPlantRenderer::WaterPlantRenderer(const Plant &plant) {
    this->init(plant);
}

void WaterPlantRenderer::init(const Plant &plant) {
    this->stem_renderer.init(plant);
    this->vertex_renderer.init(plant);
}
void WaterPlantRenderer::update(const Plant &plant) {
    this->stem_renderer.update(plant);
    this->vertex_renderer.update(plant);
}
void WaterPlantRenderer::render(Shader *shader) {
    this->stem_renderer.render(shader);
    this->vertex_renderer.render(shader);
}
