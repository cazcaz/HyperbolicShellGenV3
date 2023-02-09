#include "shellNaming.h"

ShellName::ShellName()
{
}

ShellName::~ShellName()
{
}

std::string ShellName::makeName(ShellParams &parameters)
{   
    std::string name = "Rad" +std::to_string(parameters.radius) + " Ex" + std::to_string(parameters.expansions) + " ExL" + std::to_string(parameters.extensionLength) + " mS" + std::to_string(parameters.meanStiffness) + " gS" + std::to_string(parameters.gaussStiffness) + " DC" + std::to_string(parameters.desiredCurvature) + " Str" + std::to_string(parameters.strainCoeff) + " Bs" + std::to_string(parameters.bendingStiffness);
    return name;
}
