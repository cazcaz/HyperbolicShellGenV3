#pragma once

#include <string>
#include "shellParams.h"
#include <sstream>
#include <iomanip>

class ShellName
{
public:
    ShellName();
    ~ShellName();
    std::string makeName(ShellParams &parameters);
    std::string doubleConverter(double value);

private:
};