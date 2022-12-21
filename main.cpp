#include "batchGen.h"

int main(int, char**) {
    
    //To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    std::vector<ShellParams> parameterList;
    BatchGen massCalcer;
   //test
    //Push parameters to the parameterList

    massCalcer.calculateAll(parameterList);
    return 0;
}
