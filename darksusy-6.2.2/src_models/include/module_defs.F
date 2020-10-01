ccc add here a unique integer identifier for each implemented particle model
ccc this is needed to easily switch modules at compile time, while retaining
ccc the same source code (for an example, see makefile for dsmain.F)

ccc The part below is automatically created by the script scr/update_mainmf.pl
ccc which is run from scr/preconfig.pl. Never test directly for the actual
ccc numbers below, always use the pre-compiler variables instead. The actual
ccc numbers can change when modules are added/deleted/renamed.

ccc Don't delete this comment. [MODSTART]
#define MODULE_common 1
#define MODULE_empty 2
#define MODULE_generic_decayingDM 3
#define MODULE_generic_fimp 4
#define MODULE_generic_wimp 5
#define MODULE_mssm 6
#define MODULE_silveira_zee 7
#define MODULE_vdSIDM 8
ccc Don't delete this comment. [MODEND]
