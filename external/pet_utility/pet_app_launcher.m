
function pet_app_launcher(iSubject, petFullFile)
%PET_APP
% Check if pet_utility is already installed;
ProtocolInfo = bst_get('ProtocolInfo');

% Check if pet_app is installed
customApps = matlab.apputil.getInstalledAppInfo;
if ~isempty(customApps) && strcmpi(customApps.name, 'pet_app')
    bst_progress('start', 'PET Launcher', 'Launching PET Processing Utility..');
    appinfo.status='installed';
    pet_app(petFullFile,iSubject)
    bst_progress('stop');
else
    bst_progress('start', 'PET Installer', 'Installing PET Processing Utility..');
    appFullFile = bst_fullfile(bst_get('BrainstormHomeDir'), 'external/pet_utility/pet_app.mlappinstall') ;
    appinfo = matlab.apputil.install(appFullFile);
    if strcmpi(appinfo.status, 'installed')
        % matlab.apputil.run('pet_appAPP')
        pet_app(iSubject, petFullFile);
    end
end

