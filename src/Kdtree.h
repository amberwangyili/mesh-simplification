//
// Created by YILI WANG on 6/17/19.
//

#ifndef DECIMATE_KDTREE_H
#define DECIMATE_KDTREE_H


#include <vector>
#include "Pair.h"
#include <cmath>

using namespace std;

struct  Kdnode{
    int axis;
    double split;
    std::vector<Vertex> children;
    Kdnode *lson, *rson;
    Kdnode(void):axis(-1),split(0),lson(NULL),rson(NULL){}
    inline bool is_leaf(void) const {return lson == NULL && rson == NULL;}
};

struct KdCompare{
    KdCompare(int axis = 0):axis(axis){}
    bool operator ()(const Vertex &lhs, const Vertex &rhs){
        return lhs.pos[axis]<rhs.pos[axis];
    }
    int axis;
};

class Kdtree{
public:
    typedef std::vector<Vertex>::iterator iterator;
    std::vector<Vertex> list_vertex;
    Kdnode *root;
    int max_leaf_size;

    Kdtree(const std::vector<Vertex> &list_vertex,int max_leaf_size)
            :list_vertex(list_vertex),root(NULL),max_leaf_size(max_leaf_size){
        initialize();
    }

    inline void initialize(){build(root,list_vertex.begin(),list_vertex.end(),0);}

    void build(Kdnode *&root,iterator begin, iterator end, int current){
        if(begin==end){
            return;
        }
        root = new Kdnode();

        size_t n = end-begin;
        if(n<=max_leaf_size){
            root->children = std::vector<Vertex>(begin,end);
        } else{
            std::nth_element(begin,begin+n/2,end,KdCompare(current));
            root->axis = current;
            root->split = (*(begin+n/2)).pos[current];

            int next = (current+1)%3;
            build(root->lson,begin,begin+n/2,next);
            build(root->rson,begin+n/2,end,next);
        }
    };

    void traverse(std::vector<Vertex> &res, const Kdnode *root, const Vertex &center, const double r2){
        if (root->is_leaf()){
            for(auto v = root->children.begin();v!= root->children.end();v++){
                double dist = ((*v).pos-(center).pos).dot((*v).pos-(center).pos);
                if (dist<r2) res.push_back(*v);
            }
        }else{
            bool in_left = center.pos[root->axis]<root->split;
            bool ignore = (center.pos[root->axis]-root->split)*(center.pos[root->axis]-root->split)>r2;

            if(in_left){
                traverse(res,root->lson,center,r2);
                if(!ignore) traverse(res,root->rson,center,r2);
            }else{
                traverse(res,root->rson,center,r2);
                if(!ignore) traverse(res,root->lson,center,r2);
            }
        }
    }
    inline std::vector<Vertex> find_r(const Vertex &center, const double r){
        std::vector<Vertex> res;
        double r2 = r*r;
        traverse(res,root,center,r*r);
        return res;
    }

    inline std::vector<Vertex> find_r_brute(const Vertex &center, const double r){
        std::vector<Vertex> res;
        double r2 = r*r;
        for(auto v = list_vertex.begin();v!= list_vertex.end();v++){
            double dist = ((*v).pos-(center).pos).dot((*v).pos-(center).pos);
            if (dist<r2) res.push_back(*v);
        }
        return res;

    }

};

#endif //DECIMATE_KDTREE_H
