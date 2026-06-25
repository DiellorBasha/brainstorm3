# Experimental / retired

Code retired from the toolbox, kept for reference and rebuild — NOT on the production path.

## Vortex detection + tracking (retired 2026-06-25)

`bst_vortex_persistence.m`, `bst_vortex_track.m`, `bst_vortex_link_step.m`,
`process_vortex_track.m` — the persistence-based vortex-core detector + cross-frame tracker.
Retired during the differential-Helmholtz `Compute` refactor: detection is being rebuilt from
scratch, better integrated into the new architecture (consuming `process_helmholtz('Compute')`'s
`Div`/`Curl` fields). `process_vortex_track` still references the deleted `bst_helmholtz` and is
non-functional as-is — it is a reference artifact for the rebuild, not runnable code.
