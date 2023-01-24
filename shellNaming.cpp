#include "shellNaming.h"

ShellName::ShellName()
{
}

ShellName::~ShellName()
{
}

std::string ShellName::makeName(ShellParams &parameters)
{   
    std::string name = "Rad" +std::to_string(parameters.radius) + " Ex" + std::to_string(parameters.expansions) + " ExL" + std::to_string(parameters.extensionLength) + " S" + std::to_string(parameters.stiffnessRatio) + " DC" + std::to_string(parameters.desiredCurvature) + " Str" + std::to_string(parameters.strainCoeff);
    return name;
}
