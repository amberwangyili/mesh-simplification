#include <iostream>
#include <cstdlib>
#include <iostream>
#include "Pair.h"
#include "Kdtree.h"
#include "Mesh.h"
#include<sys/time.h>
#include<unistd.h>
using namespace std;


const double thresh =  0.01;

int main(int argc, char** argv)
{
    Mesh simplifier(argv[1],thresh);
    simplifier.simplify(atof(argv[3]));
    simplifier.writer(argv[2]);
    return 0;
}
//
//int main(int argc, char** argv)
//{
////    long i = 10000000L;
////    clock_t start, mid,finish;
////    start = clock();
//    Mesh simplifier("/Users/yiliwang/Desktop/网格简化/res/orginal_pic/dragon.obj",thresh);
////    mid = clock();
////    double duration_1 = (double)(mid - start) / CLOCKS_PER_SEC;
////    printf("init time is %f s \n",duration_1);
//    simplifier.simplify(0.01);
////    finish = clock();
////    double duration_2 = (double)(finish - mid) / CLOCKS_PER_SEC;
////    printf("simp time is %f s \n",duration_2);
//    simplifier.writer("/Users/yiliwang/Desktop/网格简化/res/orginal_pic/dragon-simp.obj");
//
//    return 0;
//}



