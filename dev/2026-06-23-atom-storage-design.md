# Atom storage step — design

**Goal:** Turn the Dynamics control panel from a viewer into an interactive *atom creator*: a
**Record** action captures the extrema of the shaped field at the cursor and writes them as atoms,
tagged with exactly the `(time, band, scale, operator, region)` the panel is set to.

## Decisions (approved)
- **Scope:** record at the **cursor** (snapshot). Step the cursor + record again to accumulate.
- **Grouping:** **one group per (band + operator)**; the `charge` sign is a per-occurrence descriptor
  (sources/sinks together; vortices/antivortices together).

## Where the shaped field comes from
The linked Helmholtz figure caches, per cursor time index `iT`, `Ht = St.Cache(iT)` — the band-passed
(`GetResultsValues`, FullSources display filter) + eigen-smoothed decomposition
(`bst_helmholtz('Frame', Op, J, false)`), carrying `Vtot/Virr/Vsol`, `Phi/Psi/Fmag`. The operator
selects the scalar — and that IS the atom `Function`:

| `st.curOp` | scalar | extrema → atoms | `Function` | signed |
|---|---|---|---|---|
| Total \|J\| | `Ht.Fmag` | amplitude peaks | `magnitude` | no (max only, charge +1) |
| Φ (Irrot) | `Ht.Phi` | sources (max,+1) / sinks (min,−1) | `potential` | yes |
| Ψ (Solen) | `Ht.Psi` | vortices (max,+1) / antivortices (min,−1) | `stream` | yes |

## The detector (built fresh — replaces the stripped Helmholtz cores)
`bst_dynamics('Extrema', field, VertConn, nPeaks, signed)` → pure surface-scalar extremum finder:
- local maximum = vertex whose value ≥ all `VertConn` neighbours; take the top-`nPeaks` by value.
- if `signed`, also the top-`nPeaks` local minima (charge −1); else maxima only (charge +1).
- returns `struct('iVertex',…, 'value',…, 'charge',…)` arrays. `nPeaks` from a panel **Peaks** field
  (default 3). (`process_source_atoms`' inline `i_local_maxima` can adopt this later.)

## The Record action (panel `OnRecord`, Atoms-menu "Record at cursor")
1. Read `St = getappdata(hFig,'HelmholtzState')`, cursor `iT = bst_memory('GetTimeVector',…,'CurrentTimeIndex')`,
   `Ht = St.Cache(iT)` (recompute via `bst_helmholtz('Frame', St.Op, Jt, false)` on cache miss).
2. `Scal/Function/signed` from `st.curOp`; `Surf = getappdata(hFig,'DynamicsSurf')` (VertConn, Vertices).
3. `ex = bst_dynamics('Extrema', Scal, Surf.VertConn, nPeaks, signed)`.
4. **Find-or-create** the group keyed by `(bandName, Function)` — label `"<band> <function>"`
   (e.g. `alpha stream`); top-level **simple** group; group-level `band/bandName`=`st.curBand`,
   `scale/scaleName`=`st.curScale`, `Function`. Colour per operator.
5. **Append** one occurrence per extremum at the cursor time: `times`(+t), `vertices`, `pos`=`Surf.Vertices(v)`,
   `hemi`, `strength`=value, `charge`=sign.
6. `i_apply` (redraw markers on the Helmholtz cortex + rebuild tree) and **auto-save** to `st.file`.

## Tree refinement
A top-level **simple** group (a Record group) becomes a stack whose **selection lists its atoms in the
right pane** — flat, time-sorted, columns `time · charge · strength · vertex` — instead of one tree leaf
per occurrence (records accumulate many). Band-window (extended) groups are unchanged.

## Out of scope (later)
- Persistence/robustness gating of extrema (currently top-N by value).
- Trajectory linking across cursor steps (atoms are points now; trains later).
- `chirality` for Ψ beyond the `charge` sign; `phase` tagging from the refphase markers.

## Deliverables
1. `bst_dynamics('Extrema', …)` pure detector + a unit check.
2. `panel_bst_dynamics`: a **Peaks** field + an Atoms-menu **"Record at cursor"** → `OnRecord`
   (find-or-create (band,Function) group, append cursor extrema, refresh + auto-save).
3. Tree: simple groups list atoms in the right pane.
4. Regression: suite stays green; new check that Record on a synthetic frame appends a (band,Function)
   group with the expected occurrences/charge.
