#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>
#include <vector>
#include "shellParams.h"

class LineGraphOutputter
{
public:
    LineGraphOutputter(std::string& outputDirectory, std::string &fileName, std::string& plotTitle,std::string &xLabel, std::string &yLabel, ShellParams &parameters);
    ~LineGraphOutputter();
    void addXValues(std::vector<double>& xValues);
    void addData(std::vector<double>& yValues, std::string& legend);
    void writeData();
private:
    std::string m_outputDirectory;
    std::string m_fileName;
    std::string m_plotTitle;
    std::string m_xLabel;
    std::string m_yLabel;
    std::vector<std::string> m_legends;
    std::vector<double> m_XValues;
    std::vector<std::vector<double>> m_data;
    ShellParams m_parameters;
};