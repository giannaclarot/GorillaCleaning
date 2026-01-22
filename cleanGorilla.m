%% cleanGorilla_singleSubject.m
% Cleans Gorilla training + gen for one participant and saves cleaned CSVs to Dropbox (Gina model_fitting_gorilla).
% - Training: 240 trials -> 5 blocks (48 each), prints block accuracies
% - Gen: computes feature-away directly from raw columns feat_1..feat_8 (already in Gorilla CSV),
%        prints accuracies by feature-away (0..8) for Gen1, Gen2, and combined
%
% Outputs in workspace:
%   T_tra, T_gen_old, T_train_clean, T_gen1, T_gen2
%   trainACC, trainBlockACC (1x5)
%   gen1ACC, gen2ACC, genACC
%   gen1ACC_byAway (9x1), gen2ACC_byAway (9x1), genACC_byAway (9x1)

clear; clc;

%% Inputs
subjID = 'A-IDDRC-008';%CHANGE FOR EACH PARTICIPANT

trainCSV = fullfile('/Users/gianna/DCNL Dropbox/Vaidya Lab/Learning Mechanisms in ASD/ADULT_GORILLA/SUBJECT_DATA_RAW', subjID, [subjID '_training.csv']);
genCSV   = fullfile('/Users/gianna/DCNL Dropbox/Vaidya Lab/Learning Mechanisms in ASD/ADULT_GORILLA/SUBJECT_DATA_RAW', subjID, [subjID '_gen.csv']);

% Manually set the task version (must match the stimulus column name in the Gorilla export)
% Typical: 'Fish' / 'Bug' / 'Butterfly'
version_gorilla = 'Fish'; %CHANGE FOR EACH PARTICIPANT

% --- Button order coding sets ---
leftOrder1  = ["ORT","MIP","ULA"];
rightOrder1 = ["FID","NAX","ZEB"];

leftOrder2  = ["FID","NAX","ZEB"];
rightOrder2 = ["ORT","MIP","ULA"];


%% Output
dropboxRoot2 = '/Users/gianna/DCNL Dropbox/Lab_Shared_Scripts/Gina';
outDir = fullfile(dropboxRoot2, 'model_fitting_gorilla', 'output2', subjID);
if ~isfolder(outDir), mkdir(outDir); end

outTrain   = fullfile(outDir, [subjID '_training_clean.csv']);
outGen1    = fullfile(outDir, [subjID '_gen1_clean.csv']);
outGen2    = fullfile(outDir, [subjID '_gen2_clean.csv']);
outSummary = fullfile(outDir, [subjID '_summary_metrics.csv']);

%% =========================
%  Read + clean TRAINING
% =========================
opts = detectImportOptions(trainCSV);
opts = setvartype(opts,'Response','string');
T_tra = readtable(trainCSV, opts);
btnOrd_tra = inferButtonOrder(T_tra, leftOrder1, leftOrder2);


% Remove non-trial rows (your original filters)
if ismember('display', T_tra.Properties.VariableNames)
    T_tra(strcmp(T_tra.display,'Instructions'),:) = [];
end
if ismember('ZoneName', T_tra.Properties.VariableNames)
    T_tra(strcmp(T_tra.ZoneName,''),:)                = [];
    T_tra(strcmp(T_tra.ZoneName,'Zone2'),:)           = [];
    T_tra(strcmp(T_tra.ZoneName,'advancementZone'),:) = [];
end


% Guardrail (as you expect)
if height(T_tra) ~= 240
    fprintf('STOP: Training rows = %d (expected 240)\n', height(T_tra));
    return
end

% Sanity checks
reqCols = {'TrialNumber','ZoneType','Response','ANSWER','ReactionTime','feat_1','feat_2','feat_3','feat_4','feat_5','feat_6','feat_7','feat_8'};
for c = 1:numel(reqCols)
    if ~ismember(reqCols{c}, T_tra.Properties.VariableNames)
        error('TRAINING missing required column: %s', reqCols{c});
    end
end
if ~ismember(version_gorilla, T_tra.Properties.VariableNames)
    error('TRAINING missing stimulus column "%s" (check version_gorilla)', version_gorilla);
end

%% --- Build TRAINING clean table (1 row per trial; pull from response_button_text rows) ---
NTRIALS_TRAIN = 240;

T_train_clean = table( ...
    strings(NTRIALS_TRAIN,1), ...  % stim
    strings(NTRIALS_TRAIN,1), ...  % correct_response
    strings(NTRIALS_TRAIN,1), ...  % subj_response
    NaN(NTRIALS_TRAIN,1),... %subj_response_num
    NaN(NTRIALS_TRAIN,1), ...      % RT
    NaN(NTRIALS_TRAIN,1), ...      % ACC
    NaN(NTRIALS_TRAIN,1), NaN(NTRIALS_TRAIN,1), NaN(NTRIALS_TRAIN,1), NaN(NTRIALS_TRAIN,1), ...
    NaN(NTRIALS_TRAIN,1), NaN(NTRIALS_TRAIN,1), NaN(NTRIALS_TRAIN,1), NaN(NTRIALS_TRAIN,1), ... % feat_1..8
    NaN(NTRIALS_TRAIN,1), ...      % feaAway
    'VariableNames', {'stim','correct_response','subj_response','subj_response_num','RT','ACC', ...
                      'feat_1','feat_2','feat_3','feat_4','feat_5','feat_6','feat_7','feat_8', ...
                      'feaAway'} );

for i = 1:NTRIALS_TRAIN
    trialT = T_tra(ismember(T_tra.TrialNumber, i), :);
    if isempty(trialT), continue; end

    % remove timelimit rows
    trialT(strcmp(trialT.ZoneType,'timelimit_screen'),:) = [];

    % response rows only
    respRows = trialT(strcmp(trialT.ZoneType,'response_button_text'), :);
    if isempty(respRows), continue; end

    % Stimulus + correct response from response row
    T_train_clean.stim(i) = string(respRows.(version_gorilla){1});
    T_train_clean.correct_response(i) = string(respRows.ANSWER{1});

    % Carry features (take first response row)
    T_train_clean.feat_1(i) = respRows.feat_1(1);
    T_train_clean.feat_2(i) = respRows.feat_2(1);
    T_train_clean.feat_3(i) = respRows.feat_3(1);
    T_train_clean.feat_4(i) = respRows.feat_4(1);
    T_train_clean.feat_5(i) = respRows.feat_5(1);
    T_train_clean.feat_6(i) = respRows.feat_6(1);
    T_train_clean.feat_7(i) = respRows.feat_7(1);
    T_train_clean.feat_8(i) = respRows.feat_8(1);

    % Feature-away = # of zeros across feat_1..feat_8
    T_train_clean.feaAway(i) = sum([T_train_clean.feat_1(i),T_train_clean.feat_2(i),T_train_clean.feat_3(i),T_train_clean.feat_4(i), ...
                                    T_train_clean.feat_5(i),T_train_clean.feat_6(i),T_train_clean.feat_7(i),T_train_clean.feat_8(i)] == 0);

    % Enforce single consistent response among response rows (ignore blanks)
    uResp = unique(strtrim(string(respRows.Response)));
    uResp = uResp(uResp ~= "");

    if numel(uResp) == 1
        T_train_clean.subj_response(i) = uResp(1);
        T_train_clean.subj_response_num(i) = codeResponse(uResp(1), btnOrd_tra, ...
    leftOrder1, rightOrder1, leftOrder2, rightOrder2);


        % RT: first non-NaN among response rows
        rtVals = respRows.ReactionTime;
        rtVals = rtVals(~isnan(rtVals));
        if ~isempty(rtVals)
            T_train_clean.RT(i) = rtVals(1);
        end

        % ACC
        T_train_clean.ACC(i) = double(strcmp(T_train_clean.subj_response(i), T_train_clean.correct_response(i)));
    end
end




%% =========================
%  Read + clean GEN (split into Gen1/Gen2)
% =========================
opts = detectImportOptions(genCSV);
opts = setvartype(opts,'Response','string');
T_gen_old = readtable(genCSV, opts);
btnOrd_gen = inferButtonOrder(T_gen_old, leftOrder1, leftOrder2);
if isnan(btnOrd_gen), btnOrd_gen = btnOrd_tra; end
% Remove non-trial rows/screens
T_gen_old(cellfun(@(x) ~isempty(x), regexp(T_gen_old.ScreenName,'Instructions','once')),:)     = [];
T_gen_old(cellfun(@(x) ~isempty(x), regexp(T_gen_old.ScreenName,'FixationScreen','once')),:)  = [];
T_gen_old(strcmp(T_gen_old.ScreenName,''),:)                                                 = [];
T_gen_old(strcmp(T_gen_old.ZoneType,'continue_button'),:)                                    = [];

% Sanity checks
for c = 1:numel(reqCols)
    if ~ismember(reqCols{c}, T_gen_old.Properties.VariableNames)
        error('GEN missing required column: %s', reqCols{c});
    end
end
if ~ismember(version_gorilla, T_gen_old.Properties.VariableNames)
    error('GEN missing stimulus column "%s" (check version_gorilla)', version_gorilla);
end

% Find break screen (splits Gen1 vs Gen2)
serpInd = find(strcmp(T_gen_old.ScreenName,'BreakScreen'), 1, 'first');
if isempty(serpInd)
    error('BreakScreen not found in gen file.');
end

%% --- Build Gen1 and Gen2 tables (34 trials each; response_button_text aligned) ---
NTRIALS_GEN = 34;

for j = 1:2
    if j == 1
        old = T_gen_old(1:serpInd-1,:);
    else
        old = T_gen_old(serpInd+1:end,:);
    end

    new = table( ...
        strings(NTRIALS_GEN,1), ...
        strings(NTRIALS_GEN,1), ...
        strings(NTRIALS_GEN,1), ...
        NaN(NTRIALS_GEN,1), ...
        NaN(NTRIALS_GEN,1), ...
        NaN(NTRIALS_GEN,1), NaN(NTRIALS_GEN,1), NaN(NTRIALS_GEN,1), NaN(NTRIALS_GEN,1), ...
        NaN(NTRIALS_GEN,1), NaN(NTRIALS_GEN,1), NaN(NTRIALS_GEN,1), NaN(NTRIALS_GEN,1), ...
        NaN(NTRIALS_GEN,1), ...
        'VariableNames', {'stim','correct_response','subj_response','RT','ACC', ...
                          'feat_1','feat_2','feat_3','feat_4','feat_5','feat_6','feat_7','feat_8', ...
                          'feaAway'} );

    for i = 1:NTRIALS_GEN
        trialT = old(ismember(old.TrialNumber, i), :);
        if isempty(trialT), continue; end

        % remove timelimit rows
        trialT(strcmp(trialT.ZoneType,'timelimit_screen'),:) = [];

        % response rows only
        respRows = trialT(strcmp(trialT.ZoneType,'response_button_text'), :);
        if isempty(respRows), continue; end

        % stim + correct
        new.stim(i) = string(respRows.(version_gorilla){1});
        new.correct_response(i) = string(respRows.ANSWER{1});

        % carry features
        new.feat_1(i) = respRows.feat_1(1);
        new.feat_2(i) = respRows.feat_2(1);
        new.feat_3(i) = respRows.feat_3(1);
        new.feat_4(i) = respRows.feat_4(1);
        new.feat_5(i) = respRows.feat_5(1);
        new.feat_6(i) = respRows.feat_6(1);
        new.feat_7(i) = respRows.feat_7(1);
        new.feat_8(i) = respRows.feat_8(1);

        % feature-away
        new.feaAway(i) = sum([new.feat_1(i),new.feat_2(i),new.feat_3(i),new.feat_4(i), ...
                              new.feat_5(i),new.feat_6(i),new.feat_7(i),new.feat_8(i)] == 0);

        % response
        uResp = unique(strtrim(string(respRows.Response)));
        uResp = uResp(uResp ~= "");

        if numel(uResp) == 1
            new.subj_response(i) = uResp(1);
            new.subj_response_num(i) = codeResponse(uResp(1), btnOrd_gen, ...
    leftOrder1, rightOrder1, leftOrder2, rightOrder2);


            rtVals = respRows.ReactionTime;
            rtVals = rtVals(~isnan(rtVals));
            if ~isempty(rtVals)
                new.RT(i) = rtVals(1);
            end

            new.ACC(i) = double(strcmp(new.subj_response(i), new.correct_response(i)));
        end
    end

    if j == 1
        T_gen1 = new;
    else
        T_gen2 = new;
    end
end

%% =========================
%  ACCURACY METRICS 
% =========================

% accmean 
accmean = @(ACC) local_accmean(ACC, 'exM'); 

% validity: has response and RT not zero
validTra = T_train_clean.subj_response ~= "" & ~ismissing(T_train_clean.subj_response) & ~isnan(T_train_clean.RT) & T_train_clean.RT ~= 0;
vACC_tra = T_train_clean.ACC; vACC_tra(~validTra) = NaN;

validGen1 = T_gen1.subj_response ~= "" & ~ismissing(T_gen1.subj_response) & ~isnan(T_gen1.RT) & T_gen1.RT ~= 0;
vACC_gen1 = T_gen1.ACC; vACC_gen1(~validGen1) = NaN;

validGen2 = T_gen2.subj_response ~= "" & ~ismissing(T_gen2.subj_response) & ~isnan(T_gen2.RT) & T_gen2.RT ~= 0;
vACC_gen2 = T_gen2.ACC; vACC_gen2(~validGen2) = NaN;

% overall
trainACC = accmean(vACC_tra) * 100;
gen1ACC  = accmean(vACC_gen1) * 100;
gen2ACC  = accmean(vACC_gen2) * 100;

% training blocks (5 blocks of 48)
blockNum = 5;
triPerBlock = 240 / blockNum; % 48
trainBlockACC = cellfun(@(x) accmean(x), mat2cell(vACC_tra, ones(1,blockNum)*triPerBlock, 1))' * 100;

% gen by feature-away (0..8)
gen1ACC_byAway = NaN(9,1);
gen2ACC_byAway = NaN(9,1);
genACC_byAway  = NaN(9,1);

vACC_gen = [vACC_gen1; vACC_gen2];
feaAway_gen = [T_gen1.feaAway; T_gen2.feaAway];

for k = 0:8
    gen1ACC_byAway(k+1) = accmean(vACC_gen1(T_gen1.feaAway==k)) * 100;
    gen2ACC_byAway(k+1) = accmean(vACC_gen2(T_gen2.feaAway==k)) * 100;
    genACC_byAway(k+1)  = accmean(vACC_gen(feaAway_gen==k)) * 100;
end

%% =========================
%  Save cleaned CSVs
% =========================
writetable(T_train_clean, outTrain, 'FileType', 'text');
writetable(T_gen1, outGen1, 'FileType', 'text');
writetable(T_gen2, outGen2, 'FileType', 'text');

%% =========================
%  Save summary metrics CSV
% =========================
summary = table;
summary.subjID = string(subjID);

summary.train_overall = trainACC;
summary.train_block1  = trainBlockACC(1);
summary.train_block2  = trainBlockACC(2);
summary.train_block3  = trainBlockACC(3);
summary.train_block4  = trainBlockACC(4);
summary.train_block5  = trainBlockACC(5);

summary.gen1_overall = gen1ACC;
summary.gen2_overall = gen2ACC;

% gen overall by feature-away (0..8)
summary.gen_0feat = genACC_byAway(1);
summary.gen_1feat = genACC_byAway(2);
summary.gen_2feat = genACC_byAway(3);
summary.gen_3feat = genACC_byAway(4);
summary.gen_5feat = genACC_byAway(6);
summary.gen_6feat = genACC_byAway(7);
summary.gen_7feat = genACC_byAway(8);
summary.gen_8feat = genACC_byAway(9);

writetable(summary, outSummary, 'FileType', 'text');

%% =========================
%  Print
% =========================
fprintf('\n%s ACCURACY\n', subjID);
fprintf('Training overall: %.1f%%\n', trainACC);
fprintf('Training blocks (1..5): %.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%\n', trainBlockACC);

fprintf('\nGen1 overall: %.1f%% | Gen2 overall: %.1f%%\n', gen1ACC, gen2ACC);
fprintf('Gen (combined) ACC by feature-away 0..8:\n');
for k = 0:3
    fprintf('  %dfeat: %.1f%%\n', k, genACC_byAway(k+1));
end
for k = 5:8
    fprintf('  %dfeat: %.1f%%\n', k, genACC_byAway(k+1));
end

fprintf('\nSaved cleaned CSVs:\n%s\n%s\n%s\n', outTrain, outGen1, outGen2);
fprintf('Saved summary metrics:\n%s\n', outSummary);


%% ========= local function =========
function out = local_accmean(ACC, ACC_ind)
    if strcmp(ACC_ind,'inM')
        out = sum(ACC(~isnan(ACC))) / length(ACC);
    else % exM
        out = sum(ACC(~isnan(ACC))) / nnz(~isnan(ACC));
    end
end

function btnOrd = inferButtonOrder(T, leftOrder1, leftOrder2)
% Returns 1 or 2 based on Manipulation: LabelLeft; NaN if can't infer
btnOrd = NaN;

col = 'Manipulation: LabelLeft';
if ismember(col, T.Properties.VariableNames)
    leftLab = upper(strtrim(string(T.(col))));
    for i = 1:numel(leftLab)
        lab = leftLab(i);
        if lab == "" || ismissing(lab), continue; end
        if any(strcmp(lab, leftOrder1))
            btnOrd = 1; return
        elseif any(strcmp(lab, leftOrder2))
            btnOrd = 2; return
        end
    end
end
end

function code = codeResponse(resp, btnOrd, leftOrder1, rightOrder1, leftOrder2, rightOrder2)
% Matches your other script behavior:
% order 1: LEFT labels -> 2, RIGHT labels -> 1
% order 2: LEFT labels -> 2, RIGHT labels -> 1
code = NaN;

r = upper(strtrim(string(resp)));
if isnan(btnOrd) || r == "" || ismissing(r), return; end

if btnOrd == 1
    if any(strcmp(r, leftOrder1)),  code = 2;
    elseif any(strcmp(r, rightOrder1)), code = 1;
    end
elseif btnOrd == 2
    if any(strcmp(r, leftOrder2)),  code = 2;
    elseif any(strcmp(r, rightOrder2)), code = 1;
    end
end
end
