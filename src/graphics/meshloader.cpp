#include "graphics/meshloader.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "util/tiny_obj_loader.h"

#include <iostream>

#include <QString>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QRegularExpression>

using namespace Eigen;

bool MeshLoader::loadTriMesh(const std::string &filePath, std::vector<Vector3d> &vertices, std::vector<Vector3i> &faces)
{
    using namespace std;
    tinyobj::attrib_t attrib;
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;

    QFileInfo info(QString(filePath.c_str()));
    string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
                                info.absoluteFilePath().toStdString().c_str(), (info.absolutePath().toStdString() + "/").c_str(), true);
    if (!err.empty()) {
        cerr << err << endl;
    }

    if (!ret) {
        cerr << "Failed to load/parse .obj file" << endl;
        return false;
    }

    for (size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            Vector3i face;
            for (size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                face[v] = idx.vertex_index;

            }
            faces.push_back(face);

            index_offset += fv;
        }
    }

    for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
        vertices.emplace_back(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2]);
    }

    cout << "Loaded " << faces.size() << " faces and " << vertices.size() << " vertices" << endl;

    return true;
}

bool MeshLoader::loadTetMesh(const std::string &filepath, std::vector<Eigen::Vector3d> &vertices, std::vector<Eigen::Vector4i> &tets)
{
    QString qpath = QString::fromStdString(filepath);
    QFile file(qpath);

    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        std::cout << "Error opening file: " << filepath << std::endl;
        return false;
    }
    QTextStream in(&file);

    QRegularExpression vrxp("v (-?\\d*\\.?\\d+) +(-?\\d*\\.?\\d+) +(-?\\d*\\.?\\d+)");
    QRegularExpression trxp("t (\\d+) +(\\d+) +(\\d+) +(\\d+)");

    while(!in.atEnd()) {
        QString line = in.readLine();
        auto match = vrxp.match(line);
        if(match.hasMatch()) {
            vertices.emplace_back(match.captured(1).toDouble(),
                                  match.captured(2).toDouble(),
                                  match.captured(3).toDouble());
            continue;
        }
        match = trxp.match(line);
        if(match.hasMatch()) {
            tets.emplace_back(match.captured(1).toInt(),
                              match.captured(2).toInt(),
                              match.captured(3).toInt(),
                              match.captured(4).toInt());
        }
    }
    file.close();
    return true;
}

MeshLoader::MeshLoader()
{

}
