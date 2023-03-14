#include "shellNaming.h"

ShellName::ShellName()
{
}

ShellName::~ShellName()
{
}

std::string ShellName::makeName(ShellParams &parameters)
{   
    std::string name = "I" + std::to_string(parameters.surfaceIndex) + "Rad" +doubleConverter(parameters.radius) + "ExL" + doubleConverter(parameters.extensionLength) + "DC" + doubleConverter(parameters.desiredCurvature) + "Ls" + doubleConverter(parameters.lengthStiffness) + "Bs" + doubleConverter(parameters.bendingStiffness) + "P"+doubleConverter(parameters.period);
    return name;
}

std::string ShellName::doubleConverter(double value)
{
    std::stringstream stream;
    stream << std::scientific << std::setprecision(0) << value;
    return stream.str();
}
