%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 30/01/2020
% This script analyses the Probability Reversal task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
clearvars; clc; close all;
%% read results file
task_name = 'reversal';
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

%% define learning criterion
crit_target = 8;    % this is the target for consecutive correct responses within a single phase. This is defined based on the Murphy et al 2003 paper as the minimum number of trials that contain at least one trap trial (for 80/20 probabilities!)
txt = strcat({'Learning criterion is set to '},num2str(crit_target),{'. Do you want to proceed with this (y/n)?'});
rsp = input(txt{1},'s');
if rsp == 'n'
    rsp = input('Please enter the learning criterion [in number of trials]: ','s');
    crit_target = str2num(rsp);
end

%% loop through the items/trials of the task and save then in a struct
processed_data = struct;
phase_idx = 0;
for tr = 1:length(data.results.items)
    % find the variable 'task' that indicates when the task started and ended
    if isfield(data.results.items{tr},'task')
        if isfield(data.results.items{tr},'start')
            processed_data(1).starttime = data.results.items{tr}.start;
        else
            processed_data(phase_idx).endtime = data.results.items{tr}.end;
        end
    end
    % if new 'phase' is found then save all parameters in a new raw in the struct
    if isfield(data.results.items{tr},'phase')
        if isfield(data.results.items{tr},'start')
            trial_idx = 0;
            phase_idx = phase_idx + 1;
            processed_data(phase_idx).phase = data.results.items{tr}.phase;
        end
    end
    % find the information per trial
    if isfield(data.results.items{tr},'trial')
        trial_idx = trial_idx + 1;
        if trial_idx ~= data.results.items{tr}.trial
            fprintf('\nA mismatch was found between recorded and expected trials. Please double check the data before proceeding with the analysis.\n\n');
        end
        % check if response is missed
        if isfield(data.results.items{tr},'result') && strcmp(data.results.items{tr}.result,'timeout')
            processed_data(phase_idx).trials(trial_idx).response = NaN;
            processed_data(phase_idx).trials(trial_idx).correct = false;
            processed_data(phase_idx).trials(trial_idx).RT = NaN;
            processed_data(phase_idx).trials(trial_idx).feedback = 0;
        else
            processed_data(phase_idx).trials(trial_idx).response = data.results.items{tr}.selected;
            processed_data(phase_idx).trials(trial_idx).correct = strcmp(data.results.items{tr}.response,'correct');
            processed_data(phase_idx).trials(trial_idx).RT = data.results.items{tr}.responseTime;
            processed_data(phase_idx).trials(trial_idx).feedback = data.results.items{tr}.feedback;
        end
        processed_data(phase_idx).trials(trial_idx).isTrap = data.results.items{tr}.isTrap;
        stim = strsplit(data.results.items{tr}.stimulus,{'-',','}); % break down the stimulus to the following 4 information
        processed_data(phase_idx).trials(trial_idx).stimulus1.colour = stim{1};
        processed_data(phase_idx).trials(trial_idx).stimulus1.position = stim{2};
        processed_data(phase_idx).trials(trial_idx).stimulus2.colour = stim{3};
        processed_data(phase_idx).trials(trial_idx).stimulus2.position = stim{4};
    end
end

%% for each phase, compute the following
for ph = 1:length(processed_data)
    % number of missed responses for this phase
    results.(processed_data(ph).phase).missed = sum(isnan([processed_data(ph).trials.RT]));
    % number of incorrect responses for this phase
    results.(processed_data(ph).phase).errors = sum([processed_data(ph).trials.correct]==0);
    % average response time for all responses
    results.(processed_data(ph).phase).RT = nanmean([processed_data(ph).trials.RT]);
    % average response time for correct responses
    results.(processed_data(ph).phase).correctRT = nanmean([processed_data(ph).trials([processed_data(ph).trials.correct]).RT]);
    % number of response change after negative feedback (only if current response is correct)
    resp_change = diff([processed_data(ph).trials.correct]');
    if ph < length(processed_data)      % if there are more phases, then check for response change from the last trial of this phase to the first trial of the next phase
        resp_change = [resp_change; ~strcmp(processed_data(ph+1).trials(1).response,processed_data(ph).trials(end).response)];
    else
        resp_change = [resp_change; 0]; % if there are no other phases then response change can't be checked
    end
    results.(processed_data(ph).phase).probswitch = sum([resp_change([processed_data(ph).trials.isTrap])]==-1); % look for changes from correct to incorrect only
    % number of incorrect responses following trap trials
    idx = find([processed_data(ph).trials.isTrap]==1) + 1;
    idx(idx>length(processed_data(ph).trials)) = [];
    results.(processed_data(ph).phase).proberror = sum([processed_data(ph).trials(idx).correct]==0);
    % probability matching score
    results.(processed_data(ph).phase).probmatch = results.(processed_data(ph).phase).proberror / sum([processed_data(ph).trials.isTrap]);
    % find the first time that consecutive correct responses reach the "crit_target" value
    idx = strfind(double([processed_data(ph).trials.correct]),ones(crit_target,1)');
    if isempty(idx)
        results.(processed_data(ph).phase).trials_to_crit = NaN;
        results.(processed_data(ph).phase).postcrit_errors = NaN;
        results.(processed_data(ph).phase).maintenance_failure = NaN;
    else
        results.(processed_data(ph).phase).trials_to_crit = idx(1) + crit_target - 1; % the index starts from the start of correct responses, so we want to change this to the last response
        % how many errors after reaching the criterion
        results.(processed_data(ph).phase).postcrit_errors = sum([processed_data(ph).trials(idx(1)+crit_target:end).correct]==0);
        % maintenance failure score = percent of errors for the remaining trials after criterion was reached (the paper suggests to only count this if remaining trials are 10 or more, but that's not included here)
        results.(processed_data(ph).phase).maintenance_failure = results.(processed_data(ph).phase).postcrit_errors / (length(processed_data(ph).trials)-results.(processed_data(ph).phase).trials_to_crit);
    end
end
% number of trials until the participant updated their response after reversal
results.perseverance = find([processed_data(strcmp({processed_data.phase}.','reversal')).trials.correct],1,'first');
if isempty(results.perseverance)
    results.perseverance = NaN;
end
% combined number of response changes after negative feedback for all phases
results.probswitch = 0;
for ph = 1:length(processed_data)
    results.probswitch = results.probswitch + results.(processed_data(ph).phase).probswitch;
end

%% print out the results
fprintf('=======================================\n');
time_diff = datevec(processed_data(end).endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data(1).starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
results.duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
fprintf('Duration of task = %.2f seconds\n', results.duration);
fprintf('Perseverance after rule switching = %d trials\n', results.perseverance);
fprintf('Probability of switching after negative feedback (impulsivity) = %d times\n', results.probswitch);
for ph = 1:length(processed_data)
    fprintf('----------- Phase "%s" -----------\n',processed_data(ph).phase);
    fprintf('Number of trials = %d\n', length(processed_data(ph).trials));
    fprintf('Response Time (average) = %.2f seconds\n', results.(processed_data(ph).phase).RT/1000);
    fprintf('Missed responses = %d\n', results.(processed_data(ph).phase).missed);
    fprintf('Incorrect responses = %d\n', results.(processed_data(ph).phase).errors);
    fprintf('Trials to reach learning criterion (%d consecutive correct responses) = %d\n', crit_target, results.(processed_data(ph).phase).trials_to_crit);
    fprintf('Probability matching score = %2.1f%%\n', 100*results.(processed_data(ph).phase).probmatch);
    fprintf('Maintenance failure score = %2.1f%%\n', 100*results.(processed_data(ph).phase).maintenance_failure);
    fprintf('---------------------------------------\n');
end
fprintf('=======================================\n');

%% save results in the output file
save(strcat(path,file(1:end-4),'_results.mat'),'results','processed_data');
