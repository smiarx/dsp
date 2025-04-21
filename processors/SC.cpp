#include <SC_InterfaceTable.h>

extern void LoadTapeDelay();
extern void LoadSprings();
extern void LoadFilters();

extern InterfaceTable *ft;
InterfaceTable *ft;
PluginLoad(Processors)
{
    ft = inTable;
    LoadTapeDelay();
    LoadSprings();
    LoadFilters();
}
