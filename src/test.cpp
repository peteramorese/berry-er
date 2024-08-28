#include "berry/Operations.h"
#include "HyperRectangle.h"


using namespace BRY;
int main(int argc, char** argv) {
    
    HyperRectangle<2> h1(.5, 1.2, 3);
    HyperRectangle<2> h2(.1, 1.2, 4);

    bry_int_t b_deg = 3;
    DiagDegFilter<2> f(b_deg);
    DEBUG("Matrix w/o filter: \n" << h1.transformationMatrix(b_deg));
    NEW_LINE;
    DEBUG("Matrix w/ filter: \n" << h1.transformationMatrix(b_deg, &f));
    //HyperRectangle<2> h2(.5000000000009, 1.2, 3);

    //DEBUG("h1 < h2? " << (h1 < h2));
    //DEBUG("h2 < h1? " << (h2 < h1));
    return 0;
}