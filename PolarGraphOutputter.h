#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>
#include <vector>
#include "shellParams.h"

class PolarGraphOutputter
{
public:
    PolarGraphOutputter(std::string outputDirectory, std::string fileName, std::string plotTitle, std::string colourbarLabel, std::string cmap, ShellParams parameters, bool clipData = true);
    ~PolarGraphOutputter();
    void addRValues(std::vector<double> &rValues);
    void addData(std::vector<std::vector<double>> &polarData);
    void writeData();

private:
    std::string m_outputDirectory;
    std::string m_fileName;
    std::string m_plotTitle;
    std::string m_colourbarLabel;
    std::string m_cmap;
    bool m_clipData;
    std::vector<double> m_RValues;
    std::vector<std::vector<double>> m_data;
    ShellParams m_parameters;
};