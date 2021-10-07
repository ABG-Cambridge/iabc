%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 30/10/2020
% This script analyses the Iowa Gambling task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
clearvars; clc; close all;
%% read results file
task_name = 'iowa';
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

%% define how to split the trials
split_trials = 20;    % this is the target for consecutive correct responses within a single phase. This is defined based on the Murphy et al 2003 paper as the minimum number of trials that contain at least one trap trial (for 80/20 probabilities!)
txt = strcat({'Split trials by '},num2str(split_trials),{'. Do you want to proceed with this (y/n)?'});
rsp = input(txt{1},'s');
if rsp == 'n'
    rsp = input('Please enter the new value: ','s');
    split_trials = str2num(rsp);
end

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
        processed_data.trials(trial_idx).response = str2num(data.results.items{tr}.response(end));
        processed_data.trials(trial_idx).RT = data.results.items{tr}.responseTime;
        tmp_points = str2num(data.results.items{tr}.points);
        processed_data.trials(trial_idx).points_win = tmp_points(1);
        processed_data.trials(trial_idx).points_loss = tmp_points(2);
    end
end

%% compute net score and change in net score
results = struct();
results.nonrisky_choices = [];
results.risky_choices = [];
for tr = 1:split_trials:length(processed_data.trials)
    % advantageous choices
    results.nonrisky_choices = [results.nonrisky_choices sum([processed_data.trials(tr:tr+split_trials-1).response].'== 3) + sum([processed_data.trials(tr:tr+split_trials-1).response].'== 4)];
    % disadvantageous choices
    results.risky_choices = [results.risky_choices sum([processed_data.trials(tr:tr+split_trials-1).response].'== 1) + sum([processed_data.trials(tr:tr+split_trials-1).response].'== 2)];
end

%% print out the results
fprintf('=======================================\n');
fprintf('Number of trials = %d\n', length(processed_data.trials));
time_diff = datevec(processed_data.endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data.starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
results.duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
fprintf('Duration of task = %.2f seconds\n', results.duration);
results.RT = nanmean([processed_data.trials.RT])/1000;
fprintf('Response Time (average) = %.2f seconds\n', results.RT);
results.points_net = sum([processed_data.trials.points_win].') - sum([processed_data.trials.points_loss].');
fprintf('Total points earned = %d\n', results.points_net);
results.net = sum(results.nonrisky_choices) - sum(results.risky_choices);
fprintf('Net score (Non-risky - Risky choices) = %d\n', results.net);
results.net_change = (results.nonrisky_choices(end) - results.risky_choices(end)) - (results.nonrisky_choices(1) - results.risky_choices(1));
fprintf('Change in net score = %d\n', results.net_change);
fprintf('=======================================\n');

%% save results in the output file
save(strcat(path,file(1:end-4),'_results.mat'),'results','processed_data');
