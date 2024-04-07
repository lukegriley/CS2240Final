#pragma once

#include "graphics/shape.h"
#include "system.h"
#include "Types.h"
#include "graphics/camera.h"
#include "collider.h"


class Shader;

class Simulation
{
public:
    Simulation();

    void init();

    void update(double seconds);

    void draw(Shader *shader);

    void toggleWire();

    void loadCamera(Camera *_cam) {
        cam = _cam;
    }

    void toggleRunning();
    void forwardOnce() {
        if(!m_running) {
            once = true;
            m_running = true;
        }
    }

    void push() {
        if(m_running) {
            m_system.push(1.f);
        }
    }

    void pull() {
        if(m_running) {
            m_system.push(-1.f);
        }
    }

    void goToBreakpoint() {
        m_system.br = true;
    }
    System m_system;

private:
    Shape m_shape;

    bool equivalentFace(Vector3i f1, Vector3i f2) {
        std::vector<int> indices1 = {f1[0],f1[1],f1[2]};
        std::vector<int> indices2 = {f2[0],f2[1],f2[2]};

        std::sort(indices1.begin(), indices1.end());
        std::sort(indices2.begin(), indices2.end());

        return indices1 == indices2;
    }

    void ensureUnique(std::vector<Vector3i>& faces, std::vector<Vector3i>& tet_faces) {
        for(Vector3i& tet_face : tet_faces) {
            bool pushback = true;
            std::vector<Vector3i>::iterator it = faces.begin();
            while (it != faces.end()) {
                if (equivalentFace(tet_face,*it)) {
                    faces.erase(it);
                    pushback = false;
                } else {
                    it++;
                }
            }
            if (pushback) {
                faces.push_back(tet_face);
            }
        }
    }

    Vector3i sortedFace(Vector3i f1) {
        std::vector<int> indices1 = {f1[0],f1[1],f1[2]};

        std::sort(indices1.begin(), indices1.end());

        Vector3i res(indices1[0],indices1[1],indices1[2]);
        return res;
    }

    void initSphereCollider();

    //got this from a reference i put in the relevant readings section of my readme
    struct Vector3iHash {
        size_t operator()(const Vector3i& v) const {
            size_t hashValue = std::hash<int>()(v[0]);
            hashValue ^= std::hash<int>()(v[1]) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
            hashValue ^= std::hash<int>()(v[2]) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
            return hashValue;
        }
    };
    bool once = false;
    bool breakpoint = false;
    Camera *cam;

    Shape m_ground;
    Shape m_sphere;
    void initGround();
    bool m_running = false;
    std::vector<Vector3d> m_verts;
};
