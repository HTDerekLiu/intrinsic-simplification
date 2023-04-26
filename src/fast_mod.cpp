#include "fast_mod.h"

int fast_mod(
    const int input, 
    const int ceil)
{
    return input >= ceil ? input % ceil : input;
}