function varargout = bst_dynamics( varargin )
% BST_DYNAMICS: I/O and builders for the spatiotemporal sparse marker table.
%
% A "dynamics" table (db_template('dynamicsmat')) holds an array of Atoms
% (db_template('atom')) -- spatiotemporal sparse markers that fuse the time
% fields of Events with the space fields of Scouts, plus frequency, scale, and
% an open descriptor bag. Each Atom is a REFERENCE (time x space x band x scale
% + provenance), not a copy of the data: it is sufficient to re-derive the
% underlying source field on demand.
%
% USAGE:
%    T = bst_dynamics('New', Comment)                 % empty table
%    A = bst_dynamics('NewAtom')                       % empty atom (template)
%    T = bst_dynamics('Add', T, Atom)                  % append atom(s); keeps nAtoms in sync
%    T = bst_dynamics('Load', DynamicsFile)            % load + template-fill
%    OutFile = bst_dynamics('Save', OutFile, T)        % save table to disk
%
% SEE ALSO: db_template('atom'/'dynamicsmat'), process_source_atoms
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@

eval(macro_method);
end


%% ===== NEW TABLE =====
function T = New(Comment)
    if (nargin < 1), Comment = 'dynamics'; end
    T = db_template('dynamicsmat');
    T.Comment = Comment;
end


%% ===== NEW ATOM =====
function A = NewAtom()
    A = db_template('atom');
end


%% ===== ADD ATOM(S) =====
function T = Add(T, Atom)
    if isempty(Atom)
        return;
    end
    % Normalize each incoming atom against the template (fill missing fields,
    % drop unknown ones) so the struct array stays homogeneous.
    tmpl = db_template('atom');
    for i = 1:numel(Atom)
        a = struct_copy_fields(tmpl, Atom(i), 1);
        if isempty(T.Atoms)
            T.Atoms = a;
        else
            T.Atoms(end+1) = a;
        end
    end
    T.nAtoms = numel(T.Atoms);
end


%% ===== LOAD =====
function T = Load(DynamicsFile)
    % Phase 1: load by absolute path (the 'dynamics_' type is not yet registered
    % with file_gettype/file_fullpath; that comes with the Phase 2 tree node).
    T = load(DynamicsFile);
    T = struct_copy_fields(db_template('dynamicsmat'), T, 1);
    T.nAtoms = numel(T.Atoms);
end


%% ===== SAVE =====
function OutFile = Save(OutFile, T)
    T.nAtoms = numel(T.Atoms);
    bst_save(OutFile, T, 'v7');
    % Phase 1 has no dedicated tree node, so the table is not registered in the
    % DB; it is loadable by path via bst_dynamics('Load', ...). A tree node +
    % viewer come in Phase 2 (bst-java).
end
