#include "berry/Operations.h"
#include "HyperRectangle.h"


using namespace BRY;
int main(int argc, char** argv) {
    
    HyperRectangle<2> h1(.5, 1.2, 3);
    HyperRectangle<2> h2(.4999999999999999, 1.2, 4);
    //HyperRectangle<2> h2(.5000000000009, 1.2, 3);

    DEBUG("h1 < h2? " << (h1 < h2));
    DEBUG("h2 < h1? " << (h2 < h1));
    return 0;
}