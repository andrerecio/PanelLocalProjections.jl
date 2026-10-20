* Preamble for every cross-check do-file:
*     do "$PLP_ROOT/crosscheck/stata/common.do"
* Run do-files through crosscheck/run_do.py, which sets PLP_ROOT to the repository root.
version 19
clear all
set more off
set linesize 200
* Doubles everywhere: `generate` would otherwise store floats (~1e-7 relative error).
set type double
if `"$PLP_ROOT"' == "" {
    display as error "PLP_ROOT is not set: run this do-file through crosscheck/run_do.py"
    exit 198
}
* Project-local ado tree: installs and lookups never touch the user's PLUS/PERSONAL dirs.
sysdir set PLUS     "$PLP_ROOT/crosscheck/ado/plus"
sysdir set PERSONAL "$PLP_ROOT/crosscheck/ado/personal"
