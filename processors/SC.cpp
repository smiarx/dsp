#include "SC_PlugIn.h"

void LoadTapeDelay();
void LoadSprings();
void LoadFilters();

InterfaceTable *ft;
PluginLoad(Processors)
{
    ft = inTable;
    LoadTapeDelay();
    LoadSprings();
    LoadFilters();
}
