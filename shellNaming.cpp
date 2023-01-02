#include "shellNaming.h"

ShellName::ShellName()
{
}

ShellName::~ShellName()
{
}

std::string ShellName::makeName(ShellParams &parameters)
{   
    std::string name = "Ex" + std::to_string(parameters.expansions) + " ExL" + std::to_string(parameters.extensionLength) + " MS" + std::to_string(parameters.meanStiffness) + " GS" + std::to_string(parameters.gaussStiffness) + " DC" + std::to_string(parameters.desiredCurvature);
    return name;
}
