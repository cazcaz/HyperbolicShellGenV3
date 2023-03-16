#pragma once
#include "shellGen.h"
#include "shellParams.h"
#include <thread>
#include <vector>
#include <iostream>

class BatchGen
{
public:
    BatchGen();
    ~BatchGen();

    void calculateAll(std::vector<ShellParams> parameterList);
    void calculateBatch(std::vector<ShellParams> parameterList);
    void OnThreadShellGen(ShellParams parameters);

private:
    unsigned int m_threadCount;
    int m_surfaceCount;
    int m_totalSurfaceCount;
};