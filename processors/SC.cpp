#include <SC_InterfaceTable.h>

// extern void loadTapeDelay();
// extern void loadSprings();
extern void loadFilters();

extern InterfaceTable *ft;
InterfaceTable *ft;
PluginLoad(Processors)
{
    ft = inTable;
    // loadTapeDelay();
    // loadSprings();
    loadFilters();
}
