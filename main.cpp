#include "batchGen.h"

int main(int, char**) {
    
    //To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    std::vector<ShellParams> parameterList;
    BatchGen massCalcer;
    
    //Push parameters to the parameterList
    ShellParams parameters1;
    parameters1.expansions = 100;
    parameterList.push_back(parameters1);

    massCalcer.calculateAll(parameterList);
    return 0;
}
