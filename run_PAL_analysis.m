%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 15/06/2021
% This script analyses the Paired Associates Learning task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
% clearvars; clc; close all;
% %% read results file
% task_name = 'PAL';
ignore_practice = 1; % This will skip the results for the practice phase (=1). Turn this to =0 if you want them in the output.
% % select file using GUI
% [file,path,msg] = uigetfile(strcat('*_',task_name,'*.mat'));
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
fprintf('Task state: %s\n',data.results.state);

%% loop through the items/trials of the task and save then in a struct
processed_data = struct;
phase_idx = 0;
for tr = 1:length(data.results.items)
    % if new 'phase' is found then save all parameters in a new raw in the struct
    if isfield(data.results.items{tr},'phase')
        if isfield(data.results.items{tr},'start')
            block_idx = 0;
            phase_idx = phase_idx + 1;
            processed_data(phase_idx).phase = data.results.items{tr}.phase;
        end
    end
    % add PAL level
    if isfield(data.results.items{tr},'condition')
        if isfield(data.results.items{tr},'start') % to exclude increasing block_idx at the end of the block (condition variable appears to both start and end of the block)
            block_idx = block_idx + 1;
            processed_data(phase_idx).condition(block_idx).level = data.results.items{tr}.condition;
        end
    end
    % find the variable 'block' that indicates when a condition starts and ends
    if isfield(data.results.items{tr},'block')
        if isfield(data.results.items{tr},'start')
            trial_idx = 0;
            if data.results.items{tr}.block > 0 % increase block counter for condition repeats
                block_idx = block_idx + 1;
                processed_data(phase_idx).condition(block_idx).level = processed_data(phase_idx).condition(block_idx-1).level;
            end
            processed_data(phase_idx).condition(block_idx).starttime = data.results.items{tr}.start;
        else
            processed_data(phase_idx).condition(block_idx).endtime = data.results.items{tr}.end;
        end
    end
    % add stimuli information
    if isfield(data.results.items{tr},'stimuli')
        for i = 1:length(data.results.items{tr}.stimuli)
            trial_idx = trial_idx + 1;
            processed_data(phase_idx).condition(block_idx).trials(trial_idx).location = data.results.items{tr+i}.target;
            idx = find([data.results.items{tr}.stimuli.position].' == data.results.items{tr+i}.target);
            processed_data(phase_idx).condition(block_idx).trials(trial_idx).object = data.results.items{tr}.stimuli(idx).file;
            if strcmp(data.results.items{tr+i}.selected,'timeout')
                processed_data(phase_idx).condition(block_idx).trials(trial_idx).answer = NaN;
            else
                processed_data(phase_idx).condition(block_idx).trials(trial_idx).answer = data.results.items{tr+i}.selected;
            end
            processed_data(phase_idx).condition(block_idx).trials(trial_idx).RT = data.results.items{tr+i}.responseTime;
        end
    end
end

%% remove practice from processed_data
if ignore_practice == 1
    for ph = phase_idx:-1:1
        if strcmp(processed_data(ph).phase,'practiceMandatory') || strcmp(processed_data(ph).phase,'practiceOptional')
            processed_data(ph) = [];
            phase_idx = phase_idx - 1;
        end
    end
end

%% compute output variables
results = struct();
for ph = 1:length(processed_data)
    results.(processed_data(ph).phase) = [];
    % initialise variables
    results.(processed_data(ph).phase).RT = nan(length(processed_data(ph).condition),1);
    results.(processed_data(ph).phase).missed = nan(length(processed_data(ph).condition),1);
    results.(processed_data(ph).phase).duration = nan(length(processed_data(ph).condition),1);
    results.(processed_data(ph).phase).correct = nan(length(processed_data(ph).condition),1);
    results.(processed_data(ph).phase).level = nan(length(processed_data(ph).condition),1);
    results.(processed_data(ph).phase).memory_score = 0;
    
    % start processing in a run by run basis
    for r = 1:length(processed_data(ph).condition)
        if processed_data(ph).condition(r).level ~= length(processed_data(ph).condition(r).trials)
            fprintf('Error in expected trials! Condition is skipped for phase: %s and block: %d\n',processed_data(ph).phase,r);
            continue;
        end
        results.(processed_data(ph).phase).level(r) = processed_data(ph).condition(r).level;
        time_diff = datevec(processed_data(ph).condition(r).endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data(ph).condition(r).starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
        duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
        results.(processed_data(ph).phase).duration(r) = duration;
        results.(processed_data(ph).phase).RT(r) = mean([processed_data(ph).condition(r).trials(isfinite([processed_data(ph).condition(r).trials.answer]')).RT])/1000;
        results.(processed_data(ph).phase).missed(r) = sum(isnan([processed_data(ph).condition(r).trials.answer]));
        results.(processed_data(ph).phase).correct(r) = sum([processed_data(ph).condition(r).trials.answer].' == [processed_data(ph).condition(r).trials.location].');
        % add correct responses to memory score only for the first presentation of a condition or if the same condition is part of the design (i.e. participant responded correctly to all previous items) 
        if r == 1
            results.(processed_data(ph).phase).memory_score = results.(processed_data(ph).phase).memory_score + results.(processed_data(ph).phase).correct(r);
        else
            cond1 = results.(processed_data(ph).phase).level(r) ~= results.(processed_data(ph).phase).level(r-1);
            cond2 = ~cond1 & results.(processed_data(ph).phase).correct(r-1) == results.(processed_data(ph).phase).level(r-1);
            if cond1 || cond2
                results.(processed_data(ph).phase).memory_score = results.(processed_data(ph).phase).memory_score + results.(processed_data(ph).phase).correct(r);
            end
        end
        results.(processed_data(ph).phase).errors = sum(results.(processed_data(ph).phase).level) - sum(results.(processed_data(ph).phase).correct);
        results.(processed_data(ph).phase).num_trials = length(processed_data(ph).condition);
    end
end

%% print out the results
fprintf('=======================================\n');
for ph = 1:length(processed_data)
    fprintf('------ Results for phase: "%s" ------\n',processed_data(ph).phase);
    fprintf('Trials = %d\n', results.(processed_data(ph).phase).num_trials);
    fprintf('Response Time (average) = %.2f seconds\n', nanmean(results.(processed_data(ph).phase).RT));
    fprintf('Missed responses = %d\n', sum(results.(processed_data(ph).phase).missed));
    fprintf('Errors = %d\n', results.(processed_data(ph).phase).errors);
    fprintf('Memory score = %d\n', results.(processed_data(ph).phase).memory_score);
    fprintf('---------------------------------------\n');
end
fprintf('=======================================\n');

%% save results in the output file
save(strcat(path,file(1:end-4),'_results.mat'),'results','processed_data');
