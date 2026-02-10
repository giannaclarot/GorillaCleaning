

% Set up base paths (portable across computers/usernames)
P = local_paths();

CYJ_bigpath = fullfile(P.dropboxRoot, ...
    "Lab_Shared_Scripts", "DCNL_Members_Scripts", "Gina", "model_fitting_gorilla");


% CYJ_bigpath = char(CYJ_bigpath);
addpath(genpath(fullfile(CYJ_bigpath, "code")));

% Define conditional paths based on flags
% if isTD
%     dataGroup = 'TD';
%     subjInfoFile = 'subject_info_TD_IDDRC.csv';
% else
%     dataGroup = 'ASD';
%     subjInfoFile = 'subject_info_ASD_IDDRC.csv';
% end
% 
% % Define paths
behPath = char(fullfile(CYJ_bigpath, "data", "*"));
tableInfoPath = fullfile(CYJ_bigpath, "input", "subject_info_gorilla.csv");
timeID        = datestr(now,'yyyymmdd');

% Define output paths
outpath       = fullfile(CYJ_bigpath, "output", timeID);
% regressorPath = fullfile(outpath, "regressor"); 
                % set as '', if not wanting to save regressors
regressorPath = fullfile(''); % don't need regressors for gorilla
modelFitFile  = fullfile(outpath, "modelFit.csv");

% Ensure output directories exist
if ~exist(outpath, 'dir'), mkdir(outpath); end
if ~isempty(regressorPath) && ~exist(regressorPath, 'dir'), mkdir(regressorPath); end

%% Load inputs
behF      = g_ls(behPath);
tableInfo = readtable(tableInfoPath, 'VariableNamingRule', 'preserve');

%% Modeling
permNum = 1000;

for i = 1:length(behF)
    % Get subject ID from folder name
    [~, name] = fileparts(behF{i});

    % For debugging a particular subject.
    % if ~strcmp(name, 'LM002')
    %     continue
    % end
    % 
    % if isTD
    %     subjNum = erase(name, 'IDDRC-');   
    %     subjID = ['sub-LM' subjNum];
    %     ind = find(strcmp(tableInfo.subjID, subjID)); % checks whether subject in TableInfo
    % else
    % subjNum = name;
    % subjID = ['sub-' subjNum];
    % parts = split(name{1}, '_');
    subjID = name;
    ind = find(strcmp(tableInfo.subjID, subjID));

    if isempty(ind)
        fprintf('Skipping %s: no tableInfo entry found.\n', subjID);
        continue;
    end
    
    % Proceed if tableInfo entry exists
    % sesID = tableInfo.sesID{ind};
    % subj_sesID = tableInfo.subj_sesID{ind};
    versionL = tableInfo.version{ind};
    
    % Skip if subject already in master sheet
    if exist(modelFitFile, 'file')
        done = readtable(modelFitFile);
        if ismember(subjID, done.subjID)
            fprintf('Skipping %s — already processed.\n', subjID);
            continue;
        end
    end
       
    % Run model fitting
    disp(['Fitting model subject ' subjID ' — ' behF{i}]);
    
    [tableX_fit, tableX_weights] = IDDRC_calBeh(behF{i}, tableInfoPath, permNum, regressorPath);

    % Extract clean sesID (e.g., T1 instead of ses-T1)
    % sesID_only = extractAfter(sesID, 'ses-');  % gives 'T1'
    strategy = string(tableX_fit.strategy);
    
    % Create cleaned-up final table
    ePval    = str2double(tableX_fit.ePval);
    pPval    = str2double(tableX_fit.pPval);
    difPval  = str2double(tableX_fit.difPval);

    fitRow = table( ...
        string(subjID), ...
        string(versionL), ...
        tableX_fit.eFit, ...
        ePval, ...
        tableX_fit.pFit, ...
        pPval, ...
        round(tableX_fit.dFit, 3), ...
        difPval, ...
        string(strategy), ...
        'VariableNames', {'subjID', 'version', ...
                          'eFit', 'ePval', 'pFit', 'pPval', ...
                          'dFit', 'difPval', 'strategy'});

    finalRow = [fitRow, tableX_weights];

    if exist(modelFitFile, 'file')
            writetable(finalRow, modelFitFile, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(finalRow, modelFitFile);
    end


   % Update done list in memory (so if you rerun script in same MATLAB session it still skips)
    doneIDs(end+1,1) = subjID;

    fprintf('✅ Saved %s to %s\n', subjID, modelFitFile);
    nRun = nRun + 1;
end

fprintf('\n=== DONE ===\n');
fprintf('Ran: %d | Skipped (already processed): %d | Skipped (no tableInfo): %d\n', nRun, nSkip, nNoInfo);
fprintf('Output CSV: %s\n\n', modelFitFile);