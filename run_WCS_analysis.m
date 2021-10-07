%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 30/01/2020
% This script analyses the Wisconsin Card Sorting task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
clearvars; clc; close all;
%% read results file
task_name = 'wisconsin';
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
rule_change = [];
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
        processed_data.trials(trial_idx).rule = data.results.items{tr}.rule;
        % record the trials when the rule changes
        if (trial_idx == 1) || ~strcmp(processed_data.trials(trial_idx).rule,processed_data.trials(trial_idx-1).rule)
            rule_change = [rule_change trial_idx];
        end
        stim = strsplit(data.results.items{tr}.stimulus,'-'); % break down the stimulus to its 3 features
        processed_data.trials(trial_idx).stimulus.colour = stim{1};
        processed_data.trials(trial_idx).stimulus.shape = stim{2};
        processed_data.trials(trial_idx).stimulus.number = str2double(stim{3});
        if isfield(data.results.items{tr},'timeout')
            processed_data.trials(trial_idx).response = NaN;
            processed_data.trials(trial_idx).RT = NaN;
            processed_data.trials(trial_idx).correct = 0;
        else
            resp = strsplit(data.results.items{tr}.response,'-'); % break down the response to its 3 features
            processed_data.trials(trial_idx).response.colour = resp{1};
            processed_data.trials(trial_idx).response.shape = resp{2};
            processed_data.trials(trial_idx).response.number = str2double(resp{3});
            processed_data.trials(trial_idx).RT = data.results.items{tr}.responseTime;
            processed_data.trials(trial_idx).correct = data.results.items{tr}.isCorrect;
        end
    end
end

%% compute perseverative and non-perseverative errors
results = struct();
for f = 1:length(rule_change)
    offset = 1;
    % correct field's name if the rule has been repeated
    if isfield(results,'pers_error')
        while isfield(results.pers_error,strcat(processed_data.trials(rule_change(f)).rule,'_',num2str(offset)))
            offset = offset + 1;
        end
    end
    % select trials for the current rule only
    if f < length(rule_change)
        resp = [processed_data.trials(rule_change(f):rule_change(f+1)-1).correct].';
    else
        resp = [processed_data.trials(rule_change(f):end).correct].';
    end
    resp_deriv = diff([1; resp]); % response derivative (change)
    idx1 = find(resp_deriv==-1);  % find when participant switched to wrong from correct
    idx2 = find(resp_deriv==1);   % find when participant switched to correct from wrong
    % if first response is incorrect, count it towards the perseverative error but not as an non-perseverative error 
    if resp(1) == 0
        results.pers_error.(strcat(processed_data.trials(rule_change(f)).rule,'_',num2str(offset))) = idx2(1) - idx1(1);
        results.nonpers_error.(strcat(processed_data.trials(rule_change(f)).rule,'_',num2str(offset))) = sum(resp(idx2(1):end) == 0);
    else
        results.nonpers_error.(strcat(processed_data.trials(rule_change(f)).rule,'_',num2str(offset))) = sum(resp == 0);
        % if first is correct, still check the second one if it is incorrect
        if resp(2) == 0
            % if it is, count it towards both perseverative and non-perseverative error
            results.pers_error.(strcat(processed_data.trials(rule_change(f)).rule,'_',num2str(offset))) = idx2(1) - idx1(1);
        else
            % if it is not, then count the subsequent incorrect responses as non-perseverative errors only
            results.pers_error.(strcat(processed_data.trials(rule_change(f)).rule,'_',num2str(offset))) = 0;
        end
    end
end

%% print out the results
fprintf('=======================================\n');
fprintf('Number of trials = %d\n', length(processed_data.trials));
time_diff = datevec(processed_data.endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data.starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
results.duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
fprintf('Duration of task = %.2f seconds\n', results.duration);
results.RT = nanmean([processed_data.trials.RT])/1000;
fprintf('Response Time (average) = %.2f seconds\n', results.RT);
fprintf('Missed responses = %d\n', sum(isnan([processed_data.trials.RT])));
fields = fieldnames(results.pers_error);
for f = 1:length(fields)
    fprintf('%s perseverative errors = %d\n', strcat(upper(fields{f}(1)),fields{f}(2:end)), results.pers_error.(fields{f}));
    fprintf('%s non-perseverative errors = %d\n', strcat(upper(fields{f}(1)),fields{f}(2:end)), results.nonpers_error.(fields{f}));
end
fprintf('=======================================\n');

%% save results in the output file
save(strcat(path,file(1:end-4),'_results.mat'),'results','processed_data');
