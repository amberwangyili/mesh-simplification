//
// Created by YILI WANG on 6/17/19.
//

#ifndef DECIMATE_MESH_H
#define DECIMATE_MESH_H

#include <iostream>
#include <iomanip>
#include <opencv2/opencv.hpp>

using namespace std;


struct Vertex{
public:
    cv::Vec3d pos; //vertex coordinate
    cv::Matx44d Qmatrix; //Addictive Kpp = Qmatrix
    std::vector<size_t> adjFaces; //adjacent faces
    bool exist; //whether been contracted
    size_t rank; //index rank for kd-tree initialization
    size_t index; // index rank for final output

};


struct Face{
public:
    cv::Matx41d p; //presentation in p vector (ax+by+cz+d = 0, a^2+b^2+c^2 = 1)
    cv::Matx44d Kpmatrix; // p*p.t()
    std::vector<size_t> ves; //index for its vertices
    bool exist; // whether been degenerated

};


struct Pair{
    std::pair<size_t,size_t> v; //index for its vertices
    double deltav; //Qmatrix1 + Qmatrix2 = delta(v)
    cv::Vec3d bestv; //solution of the minimum of delta(v)
};


struct PairComp {
    bool operator () (const Pair& pa, const Pair& pb) {
        return (pa.deltav > pb.deltav); //operator for heap
    }
};
#endif //DECIMATE_MESH_H
