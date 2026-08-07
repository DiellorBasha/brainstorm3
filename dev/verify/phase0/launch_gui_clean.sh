#!/bin/bash
# Launch the desktop Brainstorm GUI against the CLEAN worktree, isolated DB.
#
# Mirrors run_matlab.sh's isolation exactly (java user.home override so
# bst_get('UserDir') resolves to an isolated .brainstorm folder, never the
# developer's real ~/.brainstorm), but:
#   - uses `matlab -desktop -r` (not `-batch`) so the MATLAB desktop and the
#     Brainstorm GUI stay open for interactive inspection instead of exiting
#     immediately;
#   - starts Brainstorm via `brainstorm('start','local')` (the GUI-mode
#     equivalent of `brainstorm('server','local')` used by run_matlab.sh) --
#     this is Brainstorm's own supported "local database" startup path: it
#     derives the DB dir from the (overridden) user.home as
#     <UserDir>/.brainstorm/local_db with zero interactive prompts, so the
#     kept `omega-tutorial-cortical-flow` protocol loads directly without
#     a "Set database folder" dialog step.
#
# Usage: ./launch_gui_clean.sh
set -euo pipefail
WT="$HOME/workspace/research/code/brainstorm3-clean"
USERDIR="${BST_USERDIR_CLEAN:-$HOME/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean}"
mkdir -p "$USERDIR"
MATLAB="/Applications/MATLAB_R2023b.app/bin/matlab"
[ -x "$MATLAB" ] || MATLAB="$(command -v matlab)"
# No trailing `exit`/`brainstorm stop`: -r (unlike -batch) leaves the desktop
# open for interactive use, which is the point of this script.
"$MATLAB" -desktop -r "java.lang.System.setProperty('user.home','$USERDIR'); cd('$WT'); brainstorm('start','local');"
