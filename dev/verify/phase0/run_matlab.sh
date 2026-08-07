#!/bin/bash
# Run a harness .m script against the CLEAN worktree's Brainstorm, isolated DB.
# Usage: ./run_matlab.sh /abs/path/to/script.m
set -euo pipefail
WT="$HOME/workspace/research/code/brainstorm3-clean"
SCRIPT="$1"
HARNESS="$HOME/workspace/research/code/brainstorm3/dev/verify/phase0"
# Isolated Brainstorm "user home": bst_get('UserDir') resolves via
# java.lang.System.getProperty('user.home'), which does NOT follow the $HOME
# env var on macOS (Java reads it from the OS user record, not the shell env).
# Without this override, `brainstorm server` loads the developer's real
# ~/.brainstorm/brainstorm.mat, whose DbVersion (from newer dev branches, per
# project memory) is newer than this older master checkout's hardcoded
# CurrentDbVersion in bst_startup.m, triggering an interactive "newer
# version" dialog that errors out in -batch mode (no user input available).
# Redirecting user.home keeps this harness fully isolated from the real
# Brainstorm user config, not just the protocol DB directory.
USERDIR="${BST_USERDIR_CLEAN:-$HOME/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean}"
mkdir -p "$USERDIR"
# DBDIR: Brainstorm protocol database directory. Using brainstorm('server','local')
# below auto-derives this as <UserDir>/.brainstorm/local_db and sets/saves it
# via bst_set internally (see bst_startup.m ~line 610) with zero interactive
# prompts -- this is Brainstorm's own supported mechanism for headless/scripted
# startup on a fresh user dir, and avoids the ordering bug in calling
# bst_set('BrainstormDbDir',...) + db_import(...) AFTER `brainstorm server` has
# already started (which is too late: on a first run with no DB dir configured
# yet, bst_startup itself blocks on an interactive folder-picker dialog before
# we ever get a chance to call bst_set).
DBDIR="${BST_DB_CLEAN:-$USERDIR/.brainstorm/local_db}"
mkdir -p "$DBDIR"
MATLAB="/Applications/MATLAB_R2023b.app/bin/matlab"
[ -x "$MATLAB" ] || MATLAB="$(command -v matlab)"
# addpath(HARNESS): harness-local utilities (tess_repair oracle) resolve without
# shadowing anything in the worktree (the clean branch has no tess_repair).
"$MATLAB" -batch "java.lang.System.setProperty('user.home','$USERDIR'); cd('$WT'); brainstorm('server','local'); addpath('$HARNESS'); addpath('$HARNESS/../phase1'); run('$SCRIPT'); brainstorm stop; exit(0);"
