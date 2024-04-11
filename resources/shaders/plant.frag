#version 330 core
out vec4 fragColor;

// Additional information for lighting
in vec4 normal_worldSpace;
in vec4 position_worldSpace;

uniform int wire = 0;
uniform vec3 color = vec3(1.0);
uniform float alpha = 1.0;

void main() {
    if (wire == 1) {
        fragColor = vec4(0.0, 0.0, 0.0, 1);
        return;
    }
    vec4 lightPos   = vec4(-2.0, 2.0, -3.0 , 1.0);
    vec3 lightColor = vec3(1.0f, alpha, 0.0f);
    vec4 lightDir   = normalize(-lightPos + position_worldSpace);
    float c = clamp(dot(-normal_worldSpace, lightDir), 0.5, 1);

    fragColor = vec4(color.r * c * lightColor[0], color.g * c * lightColor[0], color.b * c * lightColor[0], 1);
}
