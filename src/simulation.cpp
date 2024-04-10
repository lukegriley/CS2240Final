#include "simulation.h"
#include "graphics/meshloader.h"
#include <unordered_set>

#include <iostream>

using namespace Eigen;
using namespace std;

Simulation::Simulation() {}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    std::vector<Vector3d> vertices;
    std::vector<Vector4i> tets;
    if (MeshLoader::loadTetMesh("./example-meshes/cone.mesh", vertices, tets)) {
        // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
        //    of tetrahedra which are on the exterior surface of the object. Right now, this is
        //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
        //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
        //    are surface faces...
        std::vector<Vector3i> faces;
        std::unordered_set<Vector3i> uniqueFaces;
        m_verts = vertices;

        for(Vector4i t : tets) {
            std::vector<Vector3i> tet_faces({Vector3i(t[0],t[2],t[1]),
            Vector3i(t[0],t[1],t[3]),Vector3i(t[1],t[2],t[3]),Vector3i(t[0],t[3],t[2])});
            ensureUnique(faces,tet_faces);
        }

        cout << "face size: " << faces.size() <<endl;

        faces.insert(faces.end(),uniqueFaces.begin(),uniqueFaces.end());

        m_shape.init(vertices, faces, tets);
        Matrix3d scale = Matrix3d::Identity() * 0.01;
        Eigen::Quaterniond quaternion;
        quaternion = Eigen::AngleAxisd(M_PI/2.0, Eigen::Vector3d::UnitX());

        Eigen::Matrix3d rot = quaternion.toRotationMatrix();
        m_system.load(vertices,tets,Affine3d(Eigen::Translation3d(0, 2, 0)).matrix());
        Collider plane(ColliderType::COLL_PLANE);
        plane.loadPlane(1,0.01,-5.f,5.f,-5.f,5.f);
        m_system.addCollider(plane);
        // m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0)*scale.cast<float>()*rot.cast<float>()));
        m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0).cast<float>()));
        m_system.loadModel(Affine3d(Eigen::Translation3d(0, 2, 0)).matrix());
        m_system.loadCamera(cam);

    }

    initGround();
    initSphereCollider();
}

void Simulation::update(double seconds)
{
    // STUDENTS: This method should contain all the time-stepping logic for your simulation.
    //   Specifically, the code you write here should compute new, updated vertex positions for your
    //   simulation mesh, and it should then call m_shape.setVertices to update the display with those
    //   newly-updated vertices.

    // STUDENTS: As currently written, the program will just continually compute simulation timesteps as long
    //    as the program is running (see View::tick in view.cpp) . You might want to e.g. add a hotkey for pausing
    //    the simulation, and perhaps start the simulation out in a paused state.

    // Note that the "seconds" parameter represents the amount of time that has passed since
    // the last update

    if(m_running) {
        m_system.update(seconds);
        m_shape.setVertices(m_system.getNodesPositions());
        if(once) {
            m_running = false;
        }
    }

}

void Simulation::draw(Shader *shader)
{
    // m_shape.draw(shader);
    m_ground.draw(shader);
    // m_sphere.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

void Simulation::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces);

}

void Simulation::initSphereCollider()
{
    std::vector<Vector3d> sphere_verts;
    std::vector<Vector4i> sphere_tets;
    if (MeshLoader::loadTetMesh("./example-meshes/sphere.mesh", sphere_verts, sphere_tets)) {

        std::vector<Vector3i> faces;

        for(Vector4i t : sphere_tets) {
            std::vector<Vector3i> tet_faces({Vector3i(t[0],t[2],t[1]),
                                             Vector3i(t[0],t[1],t[3]),Vector3i(t[1],t[2],t[3]),Vector3i(t[0],t[3],t[2])});
            ensureUnique(faces,tet_faces);
        }
        m_sphere.init(sphere_verts, faces,sphere_tets);
        Eigen::Translation3f translation(1,0, 0);
        Eigen::Matrix3f scaling = Eigen::Matrix3f::Identity() * 1.0;
        Eigen::Affine3f combinedTransform = translation * scaling;

        m_sphere.setModelMatrix(combinedTransform);
        Collider sphere(ColliderType::COLL_SPHERE);
        sphere.loadSphere(Vector3d(1.0,0.0,0.0),1.0);
        m_system.addCollider(sphere);
    } else {
        cerr <<"Could not load sphere collider mesh" <<endl;
    }
}

void Simulation::toggleRunning() {
    if(!m_running) {
        m_running = true;
        once = false;
    } else {
        m_running = false;
    }
}
