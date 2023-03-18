#include "LineGraphOutputter.h"

LineGraphOutputter::LineGraphOutputter(std::string& outputDirectory,
                                       std::string &fileName,
                                       std::string& plotTitle,
                                       std::string &xLabel,
                                       std::string &yLabel,
                                       ShellParams &parameters) : m_outputDirectory(outputDirectory),
                                                                  m_fileName(fileName),
                                                                  m_plotTitle(plotTitle),
                                                                  m_xLabel(xLabel),
                                                                  m_yLabel(yLabel),
                                                                  m_parameters(parameters)
{
}

LineGraphOutputter::~LineGraphOutputter() {}

void LineGraphOutputter::addXValues(std::vector<double> &xValues)
{
    m_XValues = xValues;
}

void LineGraphOutputter::addData(std::vector<double> &yValues, std::string& legend)
{

    m_data.push_back(yValues);
    m_legends.push_back(legend);
}

void LineGraphOutputter::writeData()
{
    std::ofstream linePlotFile(m_outputDirectory + "/" + m_fileName + "linePlot.txt");
    linePlotFile << std::fixed << m_plotTitle << std::endl
                 << m_xLabel << "~" << m_yLabel << std::endl;
    for (int i=0; i<m_legends.size(); i++) {
        linePlotFile << m_legends[i];
        if (i != m_legends.size()-1) {
            linePlotFile << "~";
        } else {
            linePlotFile << std::endl;
        }
    }
    linePlotFile << std::fixed << "Res: " << m_parameters.resolution << std::endl
              << "Exp: " << m_parameters.expansions << std::endl
              << "Len: " << m_parameters.extensionLength << std::endl
              << "SC: " << m_parameters.springCoeff << std::endl
              << "BS: " << m_parameters.bendingStiffness << std::endl
              << "DC:" << m_parameters.desiredCurvature << std::endl << "?";
    for (int dataPointIndex = 0; dataPointIndex < m_XValues.size(); dataPointIndex++){
        linePlotFile << m_XValues[dataPointIndex] << " ";
        for (int dataIndex = 0; dataIndex < m_data.size(); dataIndex++) {
            linePlotFile << m_data[dataIndex][dataPointIndex];
            if (dataIndex != m_data.size() - 1) {
                linePlotFile << " ";
            }
        }
        if (dataPointIndex != m_XValues.size() - 1) {
            linePlotFile << ":";
        }
    }
    linePlotFile.close();
}
