#pragma once

#include <vector>
#include "Cell.h"

// Accumulate moments and write the gathered field.
void Sampling();
void Printing(const std::vector<CELL> &GlobalField);
