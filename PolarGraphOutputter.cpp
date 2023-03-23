#include "PolarGraphOutputter.h"
#include <cmath>

PolarGraphOutputter::PolarGraphOutputter(std::string outputDirectory,
                                         std::string fileName,
                                         std::string plotTitle,
                                         std::string colourbarLabel,
                                         std::string cmap,
                                         ShellParams parameters,
                                         bool clipData) : m_outputDirectory(outputDirectory),
                                                          m_fileName(fileName),
                                                          m_plotTitle(plotTitle),
                                                          m_colourbarLabel(colourbarLabel),
                                                          m_cmap(cmap),
                                                          m_parameters(parameters),
                                                          m_clipData(clipData)
{
}

PolarGraphOutputter::~PolarGraphOutputter() {}

void PolarGraphOutputter::addRValues(std::vector<double> &rValues)
{
    m_RValues = rValues;
}

void PolarGraphOutputter::addData(std::vector<std::vector<double>> &polarData)
{
    m_data = polarData;
}

void PolarGraphOutputter::writeData()
{
    std::ofstream polarPlotFile(m_outputDirectory + "/" + m_fileName + "polarPlot.txt");
    polarPlotFile << std::fixed << m_plotTitle << std::endl
                  << m_colourbarLabel << std::endl
                  << m_cmap << std::endl
                  << m_clipData << std::endl;
    polarPlotFile << std::fixed << "Res: " << m_parameters.resolution << std::endl
                  << "Res: " << m_parameters.resolution << std::endl
                  << "Exp: " << m_parameters.expansions << std::endl
                  << "Len: " << m_parameters.extensionLength / m_parameters.resolution << std::endl
                  << "BS: " << m_parameters.bendingStiffness / (m_parameters.springCoeff * std::pow(m_parameters.radius, 3)) << std::endl
                  << "DC:" << m_parameters.desiredCurvature / std::pow(m_parameters.radius, 2) << std::endl
                  << "?";
    for (int rIndex = 0; rIndex < m_RValues.size(); rIndex++)
    {
        polarPlotFile << std::fixed << m_RValues[rIndex] << " ";
        for (int polarDataIndex = 0; polarDataIndex < m_data[rIndex].size(); polarDataIndex++)
        {
            polarPlotFile << std::fixed << m_data[rIndex][polarDataIndex];
            if (polarDataIndex != m_data[rIndex].size() - 1)
            {
                polarPlotFile << ",";
            }
        }
        if (rIndex != m_data.size() - 1)
        {
            polarPlotFile << "|";
        }
    }
    polarPlotFile.close();
}
