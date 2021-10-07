%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 26/03/2021
% This script analyses the Remote Associates task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
clearvars; clc; close all;
%% read results file
task_name = 'associates';
% select file using GUI
[file,path,msg] = uigetfile(strcat('*_',task_name,'*.mat'));
% open file for reading
if msg == 0
    fprintf('* No file was selected *\n');
    return;
end
fname = strcat(path,file);
if ischar(fname) == 1
    fname = {fname};
end
% read all input files and concatenate them
data = load(fname{1});
for f = 2:length(fname)
    tmp = load(fname{f});
    data.results.items = [data.results.items; tmp.results.items];
end
if ~strcmp(data.results.taskId,task_name)
    fprintf('Task ID does not match the expected "%s" value: %s\n',task_name,data.results.taskId);
    fprintf('Please choose the correct file and run this script again!\n');
    return;
end

% make a printout of the file info (e.g. task, app version, etc)
fprintf('**** %s results are being analysed ****\n', strcat(upper(task_name(1)),task_name(2:end)));
fprintf('User: %s\n',data.results.userref);
fprintf('App Version: %s\n',data.results.sessionInfo.appVersion);
if isstruct(data.results.items)
    fprintf('Project: %s\n',data.results.items(1,1).project);
else
    fprintf('Project: %s\n',data.results.items{1,1}.project);
end
fprintf('Task state: %s\n',data.results.state);

%% loop through the items/trials of the task and save then in a struct
processed_data = struct;
trial_idx = 0;
for tr = 1:length(data.results.items)
    % find the variable 'task' that indicates when the task started and ended
    if isfield(data.results.items{tr},'task')
        if isfield(data.results.items{tr},'start')
            processed_data.starttime = data.results.items{tr}.start;
        else
            processed_data.endtime = data.results.items{tr}.end;
        end
    end
    % find the information per trial
    if isfield(data.results.items{tr},'trial')
        trial_idx = trial_idx + 1;
        if trial_idx ~= data.results.items{tr}.trial
            fprintf('\nA mismatch was found between recorded and expected trials. Please double check the data before proceeding with the analysis.\n\n');
        end
        processed_data.trials(trial_idx).response = data.results.items{tr}.userResponse;
        processed_data.trials(trial_idx).answer = data.results.items{tr}.answer;
        processed_data.trials(trial_idx).stimulus = data.results.items{tr}.stimulus;
        % ignore lowercase vs uppercase letters for counting correct responses
        processed_data.trials(trial_idx).correct = strcmpi(processed_data.trials(trial_idx).answer,processed_data.trials(trial_idx).response);
        if ~isfield(data.results.items{tr},'timeout')
            processed_data.trials(trial_idx).RT = data.results.items{tr}.responseTime;
        else
            processed_data.trials(trial_idx).RT = NaN;
        end
    end
end

%% compute and print accuracy
results = [];
fprintf('=======================================\n');
fprintf('Number of trials = %d\n', length(processed_data.trials));
time_diff = datevec(processed_data.endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data.starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
results.duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
fprintf('Duration of task = %.2f seconds\n', results.duration);
results.RT = nanmean([processed_data.trials.RT])/1000;
fprintf('Response Time (average) = %.2f seconds\n', results.RT);
fprintf('Missed responses = %d\n', sum(isnan([processed_data.trials.RT])));
fprintf('No response given = %d\n', length(find(cellfun(@isempty,{processed_data.trials.response}.'))) - sum(isnan([processed_data.trials.RT])));
% calculate performance
results.accuracy = sum([processed_data.trials.correct]) / length(processed_data.trials);
fprintf('Accuracy = %2.2f%%\n', 100*results.accuracy);
fprintf('=======================================\n');

%% save results in the output file
save(strcat(path,file(1:end-4),'_results.mat'),'results','processed_data');
