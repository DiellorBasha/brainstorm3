function [BstPetFile, sPet, Messages] = import_pet(iSubject, PetFile, FileFormat, isInteractive, isAutoAdjust, Comment, Labels)
% IMPORT_PET: Import a PET volume file in a subject of the Brainstorm database
% 
% USAGE: [BstPetFile, sPet, Messages] = import_pet(iSubject, PetFile, FileFormat='ALL', isInteractive=0, isAutoAdjust=1, Comment=[], Labels=[])
%               BstPetFiles = import_pet(iSubject, PetFiles, ...)   % Import multiple volumes at once
%
% INPUT:
%    - iSubject  : Index of the subject where to import the PET
%                  If iSubject=0 : import PET in default subject
%    - PetFile   : Full filename of the PET to import (format is autodetected)
%                  => if not specified : file to import is asked to the user
%    - FileFormat : String, one of the file formats in in_pet
%    - isInteractive : If 1, importation will be interactive (PET is displayed after loading)
%    - isAutoAdjust  : If isInteractive=0 and isAutoAdjust=1, reslice/resample automatically without user confirmation
%    - Comment       : Comment of the output file
%    - Labels        : Labels attached to this file (cell array {Nlabels x 3}: {index, text, RGB})
% OUTPUT:
%    - BstPetFile : Full path to the new file if success, [] if error
%    - sPet       : Brainstorm PET structure
%    - Messages   : String, messages reported by this function

%% ===== PARSE INPUTS =====
if (nargin < 3) || isempty(FileFormat)
    FileFormat = 'ALL';
end
if (nargin < 4) || isempty(isInteractive)
    isInteractive = 0;
end
if (nargin < 5) || isempty(isAutoAdjust)
    isAutoAdjust = 1;
end
if (nargin < 6) || isempty(Comment)
    Comment = [];
end
if (nargin < 7) || isempty(Labels)
    Labels = [];
end
% Initialize returned variables
BstPetFile = [];
sPet = [];
Messages = [];
% Get Protocol information
ProtocolInfo     = bst_get('ProtocolInfo');
ProtocolSubjects = bst_get('ProtocolSubjects');
% Default subject
if (iSubject == 0)
    sSubject = ProtocolSubjects.DefaultSubject;
% Normal subject 
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end
% Volume type
volType = 'PET';
if ~isempty(strfind(Comment, 'Import'))
    Comment = [];
end

%% ===== SELECT PET FILE =====
% If PET file to load was not defined : open a dialog box to select it
if isempty(PetFile)    
    % Get last used directories
    LastUsedDirs = bst_get('LastUsedDirs');
    % Get last used format
    DefaultFormats = bst_get('DefaultFormats');
    if isempty(DefaultFormats.PetIn)
        DefaultFormats.PetIn = 'ALL';
    end

    % Get PET file
    [PetFile, FileFormat, FileFilter] = java_getfile( 'open', ...
        ['Import ' volType '...'], ...   % Window title
        LastUsedDirs.ImportAnat, ...      % Default directory
        'multiple', 'files_and_dirs', ... % Selection mode
        bst_get('FileFilters', 'pet'), ...
        DefaultFormats.PetIn);
    % If no file was selected: exit
    if isempty(PetFile)
        return
    end
    % Expand file selection (if inputs are folders)
    PetFile = file_expand_selection(FileFilter, PetFile);
    if isempty(PetFile)
        error(['No ' FileFormat ' file in the selected directories.']);
    end
    % Save default import directory
    LastUsedDirs.ImportAnat = bst_fileparts(PetFile{1});
    bst_set('LastUsedDirs', LastUsedDirs);
    % Save default import format
    DefaultFormats.PetIn = FileFormat;
    bst_set('DefaultFormats',  DefaultFormats);
end

%% ===== LOAD PET FILE =====
isProgress = bst_progress('isVisible');
if ~isProgress
    bst_progress('start', ['Import ', volType], ['Loading ', volType, ' file...']);
end

% Load PET
isNormalize = 0;
sPet = in_pet(PetFile{1}, FileFormat, isInteractive, isNormalize);
if isempty(sPet)
    bst_progress('stop');
    return
end
% History: File name
if iscell(PetFile)
    sPet = bst_history('add', sPet, 'import', ['Import from: ' PetFile{1}]);
else
    sPet = bst_history('add', sPet, 'import', ['Import from: ' PetFile]);
end

%% ===== SAVE PET IN BRAINSTORM FORMAT =====
% Add a Comment field in PET structure, if it does not exist yet
if ~isempty(Comment)
    sPet.Comment = Comment;
    importedBaseName = file_standardize(Comment);
else
    if ~isfield(sPet, 'Comment') || isempty(sPet.Comment)
        sPet.Comment = 'PET';
    end
    % Use filename as comment
    [tmp__, importedBaseName] = bst_fileparts(PetFile);
    importedBaseName = strrep(importedBaseName, '.nii', '');
end

% Get subject subdirectory
subjectSubDir = bst_fileparts(sSubject.FileName);
% Produce a default anatomy filename
BstPetFile = bst_fullfile(ProtocolInfo.SUBJECTS, subjectSubDir, ['subjectimage_' importedBaseName '_volpet.mat']);
% Make this filename unique
BstPetFile = file_unique(BstPetFile);
% Save new PET in Brainstorm format
sPet = out_pet_bst(sPet, BstPetFile);

%% ===== REFERENCE NEW PET IN DATABASE ======
% New anatomy structure
sSubject.Anatomy(end + 1) = db_template('Anatomy');
sSubject.Anatomy(end).FileName = file_short(BstPetFile);
sSubject.Anatomy(end).Comment  = sPet.Comment;
% Default anatomy: do not change
if isempty(sSubject.iAnatomy)
    sSubject.iAnatomy = length(sSubject.Anatomy);
end
% Default subject
if (iSubject == 0)
    ProtocolSubjects.DefaultSubject = sSubject;
% Normal subject 
else
    ProtocolSubjects.Subject(iSubject) = sSubject;
end
bst_set('ProtocolSubjects', ProtocolSubjects);
% Save first PET as permanent default
if (length(sSubject.Anatomy) == 1)
    db_surface_default(iSubject, 'Anatomy', length(sSubject.Anatomy), 0);
end

%% ===== UPDATE GUI =====
% Refresh tree
panel_protocols('UpdateNode', 'Subject', iSubject);
panel_protocols('SelectNode', [], 'subject', iSubject, -1 );
% Save database
db_save();
% Unload PET (if a PET with the same name was previously loaded)
bst_memory('UnloadPet', BstPetFile);

%% ===== PET VIEWER =====
if isInteractive
    % Display PET
    view_pet(BstPetFile);
else
    if ~isProgress
        bst_progress('stop');
    end
end
