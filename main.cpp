#include <iostream>
#include <fstream>
#include <embree4/rtcore.h>
#include <embree4/rtcore_scene.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <Eigen/Dense>
#include "Camera.h"
#include "Renderer.h"

std::vector<std::string> split(std::string& input, char delimiter) {
    std::istringstream stream(input);
    std::string field;
    std::vector<std::string> result;
    while(std::getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

void read_obj(std::string filename, Object &object, std::vector<std::string> &usemtl) {
    std::vector<std::vector<int>> face_s;
    TriMesh triMesh;
    int x = 0;
    std::vector<int> number;
    std::vector<std::string> mtl_name;

    // obj file open
    std::ifstream file(filename);
    if(!file) {
        std::cout << "Error" << std::endl;
        std::exit(1);
    }
    std::string line;
    // loop
    while(std::getline(file, line)) {
        std::vector<std::string> strvec = split(line, ' ');
        std::vector<std::string> tmp;

        // first
        switch (line[0]) {
            case 'u':
                number.push_back(object.meshes.size() - 1);
                mtl_name.push_back(strvec.at(1));
                break;
            case 'o':
                object.meshes.push_back(triMesh);
                triMesh.vertices.clear();
                triMesh.vertex_normals.clear();
                triMesh.tex_coords.clear();
                triMesh.triangles.clear();
                break;

            case 'f':
                // read
                for(int n = 0; n < strvec.size() - 1; n++) {
                    tmp.push_back(strvec.at(n + 1));
                    std::vector<std::string> str = split(strvec.at(n + 1), '/');
                    face_s.push_back({std::stoi(str[0]), std::stoi(str[1]), std::stoi(str[2])});
                }
                // add
                triMesh.triangles.push_back({face_s.at(3*x).at(0), face_s.at(3*x+1).at(0), face_s.at(3*x+2).at(0)});
                x++;
                break;

            case 'v':
                // second
                switch (line[1]) {
                    // vertex
                    case ' ':
                        // read
                        for(int n = 0; n < strvec.size() - 1; n++) {
                            tmp.push_back(strvec.at(n + 1));
                        }
                        // add
                        triMesh.vertices.push_back({std::stof(tmp[0]), -std::stof(tmp[2]), std::stof(tmp[1])});
                        break;

                        // normal
                    case 'n':
                        // read
                        for(int n = 0; n < strvec.size() - 1; n++) {
                            tmp.push_back(strvec.at(n + 1));
                        }
                        // add
                        triMesh.vertex_normals.push_back({std::stof(tmp[0]), std::stof(tmp[1]), std::stof(tmp[2])});
                        break;

                    case 't':
                        // read
                        for(int n = 0; n < strvec.size() - 1; n++) {
                            tmp.push_back(strvec.at(n + 1));
                        }
                        // add
                        triMesh.tex_coords.push_back({std::stof(tmp[0]), std::stof(tmp[1])});
                        break;

                    default:
                        break;
                }

            default:
                break;
        }
    }

    object.meshes.push_back(triMesh);
    object.meshes.erase(object.meshes.begin());

    int j = 0;
    for(int i = 0; i < object.meshes.size(); i++) {
        if(i == number.at(j)) {
            usemtl.push_back(mtl_name.at(j));
            j++;
        }
        else
            usemtl.push_back("");
    }
}

void read_mtl(std::string filename, Object &object, std::vector<std::string> usemtl) {
    // mtl file open
    std::ifstream file(filename);
    if(!file) {
        std::cout << "Error" << std::endl;
        std::exit(1);
    }

    Material material;
    std::vector<Material> mtl_list;
    std::string line;
    // loop
    while(std::getline(file, line)) {
        std::vector<std::string> strvec = split(line, ' ');

        // first
        switch (line[0]) {
            case 'n':
                material.name = strvec.at(1);
                break;
            case 'K':
                // second
                switch (line[1]) {
                    // Ka
                    case 'a':
                        // read
                        material.Ka.x() = std::stod(strvec.at(1));
                        material.Ka.y() = std::stod(strvec.at(2));
                        material.Ka.z() = std::stod(strvec.at(3));
                        break;

                        // Kd
                    case 'd':
                        // read
                        material.Kd.x() = std::stod(strvec.at(1));
                        material.Kd.y() = std::stod(strvec.at(2));
                        material.Kd.z() = std::stod(strvec.at(3));

                        break;

                        // Ks
                    case 's':
                        // read
                        material.Ks.x() = std::stod(strvec.at(1));
                        material.Ks.y() = std::stod(strvec.at(2));
                        material.Ks.z() = std::stod(strvec.at(3));
                        break;

                        // Ke
                    case 'e':
                        // read
                        material.Ke.x() = std::stod(strvec.at(1));
                        material.Ke.y() = std::stod(strvec.at(2));
                        material.Ke.z() = std::stod(strvec.at(3));

                        break;

                    default:
                        break;
                }

                // d
            case 'd':
                // read
                material.d = std::stod(strvec.at(1));
                break;

                // illum
            case 'i':
                // read
                material.illum = std::stoi(strvec.at(1));
                mtl_list.push_back(material);
                break;

            case 'N':
                // second
                switch  (line[1]) {
                    // Ns
                    case 's':
                        // read
                        material.Ns = std::stod(strvec.at(1));
                        break;

                        // Ni
                    case 'i':
                        // read
                        material.Ni = std::stod(strvec.at(1));
                        break;

                    default:
                        break;
                }

            default:
                break;
        }
    }
    for(int i = 0; i < usemtl.size(); i++) {
        for(int j = 0; j < mtl_list.size(); j++) {
            if(mtl_list.at(j).name == usemtl.at(i)) {
                object.meshes.at(i).material = mtl_list.at(j);
            }
        }
    }
}

void readFile(std::string filename) {
    std::ifstream file(filename);
    if(!file) {
        std::cout << "Error" << std::endl;
        std::exit(1);
    }

    std::string line;
    // loop
    while(std::getline(file, line)) {
        std::cout << line << std::endl;
    }
}

int main() {
    std::string filename = "ice";

    Object object;
    std::vector<std::string> usemtl;

    read_obj(filename + ".obj", object, usemtl);
    read_mtl(filename + ".mtl", object, usemtl);

    RTCDevice deviceHandle = rtcNewDevice(nullptr);
    RTCScene sceneHandle = rtcNewScene(deviceHandle);
    RTCGeometry geometryHandle = rtcNewGeometry(deviceHandle, RTC_GEOMETRY_TYPE_TRIANGLE);

    int ntri = 0;
    int nvert = 0;

    for(int i = 0; i < object.meshes.size(); i++) {
        ntri += object.meshes.at(i).triangles.size();
        nvert += object.meshes.at(i).vertices.size();
    }

    float* vertices = (float*) rtcSetNewGeometryBuffer(geometryHandle,
                                                       RTC_BUFFER_TYPE_VERTEX,
                                                       0,
                                                       RTC_FORMAT_FLOAT3,
                                                       3 * sizeof(float),
                                                       nvert);

    int i = 0;
    for(int j = 0; j < object.meshes.size(); j++) {
        for(int k = 0; k < object.meshes.at(j).vertices.size(); k++) {
            vertices[3*i] = object.meshes.at(j).vertices[k][0];
            vertices[3*i + 1] = object.meshes.at(j).vertices[k][1];
            vertices[3*i + 2] = object.meshes.at(j).vertices[k][2];
            i++;
        }
    }

    unsigned* indices = (unsigned*) rtcSetNewGeometryBuffer(geometryHandle,
                                                            RTC_BUFFER_TYPE_INDEX,
                                                            0,
                                                            RTC_FORMAT_UINT3,
                                                            3 * sizeof(unsigned),
                                                            ntri);

    i = 0;
    for(int j = 0; j < object.meshes.size(); j++) {
        for(int k = 0; k < object.meshes.at(j).triangles.size(); k++) {
            indices[3*i] = object.meshes.at(j).triangles[k][0] - 1;
            indices[3*i + 1] = object.meshes.at(j).triangles[k][1] - 1;
            indices[3*i + 2] = object.meshes.at(j).triangles[k][2] - 1;
            i++;
        }
    }

    rtcCommitGeometry(geometryHandle);
    rtcAttachGeometry(sceneHandle, geometryHandle);
    rtcReleaseGeometry(geometryHandle);
    rtcCommitScene(sceneHandle);

    const clock_t start = clock();
    const Eigen::Vector3d campos(0, 0, 80);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;
    const Camera camera(campos, camdir, 540, 4.0 / 3.0, 60, 45);

    const Renderer renderer(sceneHandle, object, camera, Color(0.1, 0.1, 0.1));

    const unsigned int samples = 1e4;

    const auto image2 = renderer.directIlluminationRender(samples).apply_gamma_correction();


    std::cout << filename << "\nsample: " << samples << std::endl;

    image2.save(filename + ".png");

    rtcReleaseScene(sceneHandle);
    rtcReleaseDevice(deviceHandle);

    const clock_t end = clock();
    std::cout << "rendering = " << (end - start) / CLOCKS_PER_SEC << " sec" << std::endl;

    return 0;
}