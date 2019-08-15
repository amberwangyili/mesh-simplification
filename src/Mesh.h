//
// Created by YILI WANG on 6/17/19.
//

#ifndef DECIMATE_SIMPLIFIER_H
#define DECIMATE_SIMPLIFIER_H


#include "Kdtree.h"
#include "Pair.h"
#include <iostream>
#include <iomanip>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

class Mesh{
public:
    vector<Face> FaceList;
    vector<Vertex> VertexList;
    priority_queue<Pair,vector<Pair>,PairComp> PairHeap;
public:

    void reader(const string& infile_name);
    void writer(const string& outfile_name);
    void init_vertex();
    void init_face();
    void init_pair(double threshold);

    void update_face_matrix(Face &face);
    void update_vertice_matrix(Vertex &vertex);
    void update_pair(Pair &pair);

    bool isBorder(size_t vertex_r);
    bool isFlipped(size_t face_r,Pair& pair);

    Mesh(const string& infile_name,double threshold = 0.01){
        reader(infile_name);
        init_face(); //initialization
        init_vertex();
        init_pair(threshold);
    }

    void simplify(double rate);

};





void Mesh::reader(const string &infile_name) {
    FILE *inobj = fopen(infile_name.c_str(), "r");
    char *buffer = new char[1024];
    int point_rank = 0;
    while (!feof(inobj)) {
        fgets(buffer, 1024, inobj);
        if (buffer[0] == 'v') {
            Vertex point;
            point.exist = true;
            sscanf(buffer, "v %lf %lf %lf", &point.pos[0], &point.pos[1], &point.pos[2]);
            point.rank = point_rank++;
            VertexList.push_back(point);
//            cout<<"point "<<point.pos[0]<<" "<<point.pos[1]<<" "<<point.pos[2]<<" "<<"rank: "<<point.rank<<endl;
        } else if (buffer[0] == 'f') {
            Face face;
            size_t fv[3];
            face.exist = true;
            sscanf(buffer, "f %lu %lu %lu", &fv[0], &fv[1], &fv[2]);
            for (int i = 0; i < 3; ++i) {
                fv[i] -=1;
                face.ves.push_back(fv[i]);
            }
            FaceList.push_back(face);
        }
    }
    fclose(inobj);
    delete  []buffer;
}



void Mesh::writer(const string& outobjfile) {
    ofstream outobj(outobjfile);
    size_t final_index = 1; //obj start from one!
    for (auto i = VertexList.begin(); i != VertexList.end(); ++i) {
        if ((*i).exist) {
            (*i).index = final_index;
            ++final_index;
            outobj << "v " << (*i).pos[0] << ' ' << (*i).pos[1] << ' ' << (*i).pos[2] << endl;
        }
    }
    for (auto i = FaceList.begin(); i != FaceList.end(); ++i) {
        if ((*i).exist) {
            vector<size_t> v((*i).ves.begin(), (*i).ves.end());
            outobj << "f " << VertexList.at(v[0]).index << ' ' << VertexList.at(v[1]).index << ' ' << VertexList.at(v[2]).index << endl;
        }
    }
    outobj.close();
}


void Mesh::init_face() {
    for(auto i = FaceList.begin(); i!= FaceList.end();i++){
        //get coordinates of a face
        vector<cv::Vec3d> vertices;
        //for a face = (a,b,c)
        for(auto j = (*i).ves.begin(); j != (*i).ves.end();++j){
            vertices.push_back(VertexList.at(*j).pos); //vertices = (xa,ya,za)
        }
        //get its kpmat
        cv::Vec3d norm = normalize((vertices.at(1)-vertices.at(0)).cross(
                vertices.at(2)-vertices.at(0)));
        (*i).p = Matx41d::zeros();
        for (int k = 0; k < 3; ++k)
            (*i).p(k) = norm[k];
        (*i).p(3) = - (norm.dot(vertices.at(0)));
        (*i).Kpmatrix = (*i).p * (*i).p.t();

    }
}

void Mesh::init_vertex() {
    for (auto i = FaceList.begin(); i != FaceList.end(); i++) {
        Face& face = (*i); //face k = (a,b,c)
        for (auto fv = face.ves.begin(); fv != face.ves.end(); fv++) {
            VertexList.at(*fv).adjFaces.push_back(i - FaceList.begin()); //a.push_back(k)
        }
    }

    for (auto i = VertexList.begin(); i != VertexList.end(); ++i) { //for vertex 1,2,3,..,n
        (*i).Qmatrix = Matx44d::zeros(); //init
        for (auto j = (*i).adjFaces.begin(); j!= (*i).adjFaces.end() ; ++j) {
            (*i).Qmatrix += FaceList.at(*j).Kpmatrix; //Q=Kp1+..+Kpk
        }
    }
}


void Mesh::init_pair(double threshold) {
    set<pair<size_t, size_t> > pairset;

    Kdtree tree = Kdtree(VertexList,30); //initialize a kdtree for all vertices

    for (auto j = VertexList.begin();  j!= VertexList.end() ; ++j) {
        vector<Vertex> res = tree.find_r(*j,threshold); //find all nearby vertices around verteix j
        for (auto k = res.begin(); k!=res.end(); ++k) {
            size_t v1 =(*j).rank, v2 = (*k).rank;
            if(v1>v2) swap(v1,v2);
            if(pairset.find(make_pair(v1,v2))==pairset.end()&&v1!=v2){ //if not been added to pair heap before (use set to justify)
                pairset.insert(make_pair(v1,v2)); //
                Pair pair;
                pair.v = make_pair(v1,v2);
                update_pair(pair); //initialize pair's delta, bestv

                PairHeap.push(pair);
            }
        }
    }
    for (auto i = FaceList.begin(); i != FaceList.end(); ++i) {

        vector<size_t> vers;  // verteices = (a,b,c) belongs to face i
        for (auto v = (*i).ves.begin(); v !=(*i).ves.end(); ++v) {
            vers.push_back(*v);
        }

        for (int ver = 0; ver < 3; ++ver) { //for edge ab bc ca
            size_t v1 = vers[ver], v2 = vers[(ver + 1) % 3];
            if (v1 > v2)
                swap(v1, v2);
            if (pairset.find(make_pair(v1, v2)) == pairset.end()&&v1!=v2) { //check whether already in the set
                pairset.insert(make_pair(v1, v2));
                Pair pair;
                pair.v = make_pair(v1, v2);
                update_pair(pair); //initialize pair's delta, bestv

                PairHeap.push(pair);
            }
        }



    }
}

void Mesh::update_face_matrix(Face &face) {

    vector<Vec3d> tempv;
    for (auto i = face.ves.begin(); i != face.ves.end(); ++i) {
        tempv.push_back(VertexList.at(*i).pos); //get coordinate of the vertices of face i = (a,b,v)
    }
    Vec3d abc = normalize((tempv.at(1) - tempv.at(0)).cross(tempv.at(2) - tempv.at(0)));
    face.p = Matx41d::zeros();
    for (int k = 0; k < 3; ++k)
        face.p(k) = abc[k];
    face.p(3) = - (abc.dot(tempv.at(0)));
    face.Kpmatrix = face.p * face.p.t();
}

void Mesh::update_vertice_matrix(Vertex &vertex) { //set Qmatrix
    vertex.Qmatrix = Matx44d::zeros();
    set<size_t> neighbor_faces;
    int count = 0;
//    for (auto i = vertex.adjFaces.begin(); i != vertex.adjFaces.end(); ++i) {
//        neighbor_faces.insert(*i);
//    }
//    for (auto i = neighbor_faces.begin(); i != neighbor_faces.end(); ++i) {
//        vertex.Qmatrix += FaceList.at(*i).Kpmatrix;
//    }
    for (auto i = vertex.adjFaces.begin(); i != vertex.adjFaces.end(); ++i) {
        vertex.Qmatrix += FaceList.at(*i).Kpmatrix;
    }
}


void Mesh::update_pair(Pair &pair) {

    Matx44d Qvmatrix = VertexList.at(pair.v.first).Qmatrix + VertexList.at(pair.v.second).Qmatrix; //Q1+Q2
    Matx44d derivative = Matx44d::zeros();

    //|q11 q12 q13 q14|
    //|q12 q22 q23 q24| * bestv = res
    //|q13 q24 q33 q34|
    //| 0   0   0   1 |
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i > j)
                derivative(i, j) = derivative(j, i);
            else
                derivative(i, j) = Qvmatrix(i, j);
        }
    }
    derivative(3, 0) = derivative(3, 1) = derivative(3, 2) = 0;
    derivative(3, 3) = 1;

    Matx41d res = Matx41d::zeros(); res(3) = 1;
    Matx41d bestv = Matx41d::zeros();

    int judge = solve(derivative, res, bestv);
    if (judge == 1) {}
    else {
        solve(derivative, res, bestv, DECOMP_SVD); //if not invertible, then compute its pseudo-inverse [ref: http://graphics.stanford.edu/courses/cs468-08-winter/slides/jerry1.pdf]
//        pair.bestv = (VertexList.at(pair.v.first).pos+VertexList.at(pair.v.second).pos)/2;
//        cout<<"not invertible"<<endl;
    }
    for(int i = 0; i < 3; ++i)
        pair.bestv[i] = bestv(i);
    pair.deltav = (bestv.t() * Qvmatrix * bestv)(0);



}



void Mesh::simplify(double rate) {

    size_t former_face_num = FaceList.size();
    size_t after_face_num = former_face_num;
    
    while ((double)after_face_num / (double)former_face_num > rate) {

        Pair contract_pair = PairHeap.top();
        PairHeap.pop();

        //check border condition
        bool border2 = isBorder(contract_pair.v.second);
        bool border1 = isBorder(contract_pair.v.first);

        if(border1&&border2){
//            cout<<"both are borders"<<endl;
            contract_pair.bestv =  (VertexList.at(contract_pair.v.first).pos+VertexList.at(contract_pair.v.second).pos)/2;
        }
        if(border1&&!border2){
//            cout<<"border1"<<endl;
            contract_pair.bestv =  (VertexList.at(contract_pair.v.first)).pos;
        }
        if(!border1&&border2){
//            cout<<"border2"<<endl;
            contract_pair.bestv =  (VertexList.at(contract_pair.v.second)).pos;
        }

        if (!VertexList.at(contract_pair.v.first).exist || !VertexList.at(contract_pair.v.second).exist) continue;



        VertexList.at(contract_pair.v.first).exist = false; //set pair vertices invalid
        VertexList.at(contract_pair.v.second).exist = false;
        
        VertexList.push_back(Vertex());
        size_t new_ver = VertexList.size() - 1; //new vertices index
        VertexList.at(new_ver).pos = contract_pair.bestv; //set its coordinate
        VertexList.at(new_ver).exist = true; //set its state

        set<size_t> neighbor_faces; //collect all the face relevant to the pair
        set<size_t> neighbor_vers; //collect all the vertices relevant to the pair

        neighbor_faces.insert(VertexList.at(contract_pair.v.first).adjFaces.begin(), VertexList.at(contract_pair.v.first).adjFaces.end());
        neighbor_faces.insert(VertexList.at(contract_pair.v.second).adjFaces.begin(), VertexList.at(contract_pair.v.second).adjFaces.end());

        for (auto face = neighbor_faces.begin(); face != neighbor_faces.end(); ++face) {
            neighbor_vers.insert(FaceList.at(*face).ves.begin(), FaceList.at(*face).ves.end());
        }

        neighbor_vers.erase(contract_pair.v.first);
        neighbor_vers.erase(contract_pair.v.second);

        //deal with degenerate face as well as all relevant faces' index 
        //update its relevant vertices adjacent face list as well as the new vertices adjacent faces
        for (auto f = neighbor_faces.begin(); f != neighbor_faces.end(); ++f) {

            //if not flipped
//            if(isFlipped(*f,contract_pair)){continue;}
            for (auto v = FaceList.at(*f).ves.begin(); v != FaceList.at(*f).ves.end(); ++v) {
                if ((*v) == contract_pair.v.first || (*v) == contract_pair.v.second)
                    (*v) = new_ver; //substitude v1 v2 for new vertices index for all relevant faces
            }
            
            set<size_t> f_vertices; //check for invalid faces (degenerate case)
            f_vertices.insert(FaceList.at(*f).ves.begin(), FaceList.at(*f).ves.end());
            if (f_vertices.size() != 3) {
                FaceList.at(*f).exist = false;
                //delete it in relevant vertices adjacent faces
                for (auto vv = FaceList.at(*f).ves.begin(); vv != FaceList.at(*f).ves.end(); ++vv) {
                    if ((*vv) != new_ver) {
                        auto target = find(VertexList.at(*vv).adjFaces.begin(), VertexList.at(*vv).adjFaces.end(), *f);
                        VertexList.at(*vv).adjFaces.erase(target);
                    }
                }
                --after_face_num;

            }
            else {
                //add this face to the adjacent face of the new vertices
                update_face_matrix(FaceList.at(*f));
                VertexList.at(new_ver).adjFaces.push_back(*f);
            }
        }
        update_vertice_matrix(VertexList.at(new_ver));
        
        //deal with the relevant vertices (after dealing with the faces!)
        vector<size_t> updated_vers; //collect new vertices
        updated_vers.push_back(new_ver); 
        for (auto vr = neighbor_vers.begin(); vr != neighbor_vers.end(); ++vr) {
            update_vertice_matrix(VertexList.at(*vr));//update its matrix
            VertexList.at(*vr).exist = false; //set its origin as invalid
            VertexList.push_back(Vertex()); 
            size_t nv = VertexList.size() - 1;
            VertexList.at(nv) = VertexList.at(*vr); //copy relevant info from its origin
            VertexList.at(nv).exist = true;
            for (auto fr = VertexList.at(*vr).adjFaces.begin(); fr != VertexList.at(*vr).adjFaces.end(); ++fr) {
                auto target = find(FaceList.at(*fr).ves.begin(), FaceList.at(*fr).ves.end(), *vr);  
                (*target) = nv; //set its relevant face index to new index
            }
            updated_vers.push_back(nv); 
        }
        
        //deal with relevant pair
        set<pair<size_t, size_t> > updated_pair; 
        for (auto v = updated_vers.begin(); v != updated_vers.end(); ++v) {
            set<size_t> neighbor_v;
            //add the edge that is relevant to the updated point to the pair
            for (auto f = VertexList.at(*v).adjFaces.begin(); f != VertexList.at(*v).adjFaces.end(); ++f) {
                neighbor_v.insert(FaceList.at(*f).ves.begin(), FaceList.at(*f).ves.end());
            }

            for (auto near_v = neighbor_v.begin(); near_v != neighbor_v.end(); ++near_v) {
                if (*near_v != *v) {
                    auto pair = make_pair(*near_v, *v);
                    if (pair.first > pair.second)
                        swap(pair.first, pair.second);
                    if (updated_pair.find(pair) == updated_pair.end()) {
                        Pair new_pair;
                        new_pair.v = pair;
                        update_pair(new_pair);


                        PairHeap.push(new_pair);


                        updated_pair.insert(pair);
                    }
                }
            }
        }
    }
}

bool Mesh::isBorder(size_t vertex_r) {

    map<size_t, int> times;
    Vertex& vertex = VertexList.at(vertex_r);
    for (auto f = vertex.adjFaces.begin(); f != vertex.adjFaces.end(); ++f) {
        for (auto v = FaceList.at(*f).ves.begin(); v != FaceList.at(*f).ves.end(); ++v) {
            if (*v != vertex_r) {
                auto target = times.find(*v);
                if (target == times.end()) {
                    times.insert(make_pair(*v, 1));
                }
                else {
                    ++(*target).second;
                }
            }
        }
    }
    for (auto m = times.begin(); m != times.end(); ++m) {
        if ((*m).second == 1)
            return true;
    }
    return false;
}
//bool Mesh::isFlipped(size_t face_r,Pair& pair) {
//
//
//    vector<Vec3d> *check;
//    int appear_times = 0;
//    for (auto v = FaceList.at(face_r).ves.begin(); v != FaceList.at(face_r).ves.end(); ++v) {
//        if ((*v) == pair.v.first || (*v) == pair.v.second){
//            check->push_back(pair.bestv);
//            appear_times++;
//        }
//        else check->push_back(VertexList.at(*v).pos);
//    }
//    if(appear_times==1){
//        Vec3d norm = normalize((check->at(1) - check->at(0)).cross(check->at(2) - check->at(0)));
//        Vec3d original_norm = Vec3d(FaceList.at(face_r).p(0),FaceList.at(face_r).p(1),FaceList.at(face_r).p(2));
//        if(original_norm.dot(norm)<0){
//            cout<<"flipped"<<endl;
//            delete check;
//            return true;
//        }
//
//    }
//    delete check;
//    return false;
//
//}

#endif //DECIMATE_SIMPLIFIER_H
