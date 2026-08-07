# Phase 0 follow-ups (deferred minors from reviews)

Clean-branch polish before any real upstream PR:
- process_segment_freesurfer.m:286 interactive flow now shows the new method dialog (scope leak) — decide desired behavior or pin reducepatch there.
- process_import_bids: nVertices unvalidated on icosphere->reducepatch fallback for non-FS subjects.
- 'snaps to ico grid' comments overstate tess_sphere's snap list (cosmetic).
- ResolveAnatDownsample DEFERRED comment style (cosmetic).

Harness:
- run_matlab.sh dead DBDIR/BST_DB_CLEAN variable; misleading macOS user.home comment.
- launch_gui_clean.sh untested (mirrors tested runner isolation).
