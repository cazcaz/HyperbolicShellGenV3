#include "shellNaming.h"

ShellName::ShellName()
{
}

ShellName::~ShellName()
{
}

std::string ShellName::makeName(ShellParams &parameters)
{   
    std::string name = "Rad" +doubleConverter(parameters.radius) + "Ex" + doubleConverter(parameters.expansions) + "ExL" + doubleConverter(parameters.extensionLength) + "mS" + doubleConverter(parameters.meanStiffness) + "gS" + doubleConverter(parameters.gaussStiffness) + "DC" + doubleConverter(parameters.desiredCurvature) + "Str" + doubleConverter(parameters.strainCoeff) + "Bs" + doubleConverter(parameters.bendingStiffness) + "P"+doubleConverter(parameters.period);
    return name;
}

std::string ShellName::doubleConverter(double value)
{
    std::stringstream stream;
    stream << std::scientific << std::setprecision(0) << value;
    return stream.str();
}
