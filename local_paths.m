function P = local_paths()
%LOCAL_PATHS  Machine-independent paths (works across different usernames).
% Data stays on Dropbox; code stays in the Git repo.

homeDir = char(java.lang.System.getProperty('user.home'));  % same as getenv('HOME')

% ---- Find your Dropbox root (try common names first) ----
candidates = {
    fullfile(homeDir, "DCNL Dropbox")
    fullfile(homeDir, "Dropbox")
    fullfile(homeDir, "Dropbox (Personal)")
    fullfile(homeDir, "Dropbox (Work)")
};

dropboxRoot = "";
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        dropboxRoot = candidates{i};
        break
    end
end

if dropboxRoot == ""
    error("Could not find Dropbox folder in your home directory. Update candidates in local_paths.m.");
end

% ---- Set your project-specific Dropbox locations here ----
P.dropboxRoot = dropboxRoot;

% Example: where your Gorilla DATA lives on Dropbox (edit this once to match your real structure)
P.gorillaDataDir = fullfile(dropboxRoot, "Lab_Shared_Scripts", "Gorilla_Data");  % <-- change me

% Example: where you want outputs written (can be inside repo, but you already ignore /output/)
P.outputDir = fullfile(fileparts(mfilename('fullpath')), "output");
if ~isfolder(P.outputDir), mkdir(P.outputDir); end
end
