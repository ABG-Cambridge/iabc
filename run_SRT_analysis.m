%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 26/04/2021
% This script analyses the Serial Response Time task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
clearvars; clc; close all;
%% read results file
task_name = 'SRT';
ignore_practice = 1; % This will skip the results for the practice phase (=1). Turn this to =0 if you want them in the output.
remove_short_runs = 100; % Set this number to the minimum acceptable trials (shorter runs will be removed!)
% select files using GUI (it's recommended to select ALL sessions of a single user to provide a more accurate model estimation)
[files,path,msg] = uigetfile(strcat('*_',task_name,'*.mat'),'MultiSelect','on');
% open file for reading
if msg == 0
    fprintf('* No file was selected *\n');
    return;
end
fname = strcat(path,files);
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
fprintf('**** %s results are being analysed ****\n',task_name);
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
    % find the level and design information
    if isfield(data.results.items{tr},'condition') % revisit this if multiple conditions are allowed per phase!!!
        processed_data(phase_idx).level = data.results.items{tr}.level;
        processed_data(phase_idx).symbolMapping = char(strsplit(data.results.items{tr}.symbolMapping,','))';
        processed_data(phase_idx).design = data.results.items{tr}.transitionFreqs;
        % convert transition frequencies to a numerical table
        tmp = strsplit(processed_data(phase_idx).design,{'[',']'});
        processed_data(phase_idx).design = [];
        switch processed_data(phase_idx).level
            % for level-0 read the design as columns
            case 0
                for i = 1:length(tmp)
                    if ~isempty(tmp{i})
                        processed_data(phase_idx).design = [processed_data(phase_idx).design str2num(tmp{i})];
                    end
                end
            % for level-1 read the design as rows
            case 1
                for i = 1:length(tmp)
                    if ~isempty(tmp{i})
                        processed_data(phase_idx).design = [processed_data(phase_idx).design; str2num(tmp{i})];
                    end
                end
            % for level-2 read the design as rows and convert it into a 16x4 matrix
            case 2
                for i = 1:length(tmp)
                    if ~isempty(tmp{i})
                        processed_data(phase_idx).design = [processed_data(phase_idx).design; str2num(tmp{i})];
                    end
                end
                processed_data(phase_idx).design = [sum(processed_data(phase_idx).design(:,1:4),2) sum(processed_data(phase_idx).design(:,5:8),2) sum(processed_data(phase_idx).design(:,9:12),2) sum(processed_data(phase_idx).design(:,13:16),2)];
        end
    end
    % find the variable 'block' that indicates when a run starts and ends
    if isfield(data.results.items{tr},'block')
        if isfield(data.results.items{tr},'start')
            trial_idx = 0;
            block_idx = block_idx + 1;
            if isfield(data.results.items{tr},'sequence')
                processed_data(phase_idx).runs(block_idx).sequence = data.results.items{tr}.sequence' + 1;
            else
                processed_data(phase_idx).runs(block_idx).sequence = data.results.items{tr}.sequenceItems' + 1;
            end
            processed_data(phase_idx).runs(block_idx).starttime = data.results.items{tr}.start;
        else
            processed_data(phase_idx).runs(block_idx).endtime = data.results.items{tr}.end;
            processed_data(phase_idx).runs(block_idx).avg_RT = data.results.items{tr}.blockPerformanceIndex;
        end
    end
    % find the information per trial (if PI and distributions above are correct, then you might not need this!!!)
    if isfield(data.results.items{tr},'trial')
        trial_idx = trial_idx + 1;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).correct = strcmp(data.results.items{tr}.result,'correct');
        if isfield(data.results.items{tr},'response')
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).response = find(data.results.items{tr}.response == processed_data(phase_idx).symbolMapping);
        else
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).response = NaN;
        end
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).RT = data.results.items{tr}.responseTime;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).tJitter = data.results.items{tr}.jitter;
        % classify trials as "High" probable, "Low" probable or "None" (this applies to the first trials only)
        if processed_data(phase_idx).level == 1
            if trial_idx == 1
                processed_data(phase_idx).runs(block_idx).trials(trial_idx).prob = 'N';
            else
                switch processed_data(phase_idx).design(processed_data(phase_idx).runs(block_idx).sequence(trial_idx-1),processed_data(phase_idx).runs(block_idx).sequence(trial_idx))
                    case max(processed_data(phase_idx).design(processed_data(phase_idx).runs(block_idx).sequence(trial_idx-1),:))
                        processed_data(phase_idx).runs(block_idx).trials(trial_idx).prob = 'H';
                    case min(nonzeros(processed_data(phase_idx).design(processed_data(phase_idx).runs(block_idx).sequence(trial_idx-1),:)))
                        processed_data(phase_idx).runs(block_idx).trials(trial_idx).prob = 'L';
                end
            end
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

%% remove short runs
for ph = phase_idx:-1:1
    if ~strcmp(processed_data(ph).phase,'test')
        continue;
    end
    for bl = length(processed_data(ph).runs):-1:1
        if length(processed_data(ph).runs(bl).trials) < remove_short_runs
            processed_data(ph).runs(bl) = [];
        end
    end
    % if a phase becomes empty then remove
    if length(processed_data(ph).runs) == 0
        processed_data(ph) = [];
    end
end

%% combine data from multiple sessions if phase is the same
phases = unique({processed_data.phase});
phases_to_remove = [];
for ph = 1:length(phases)
    % find all results from the same phase to combine them
    idx = find(strcmp({processed_data.phase},phases{ph})==1);
    ref_struct = processed_data(idx(1));
    ref_struct.runs = []; % remove the run details so the two structs can be compared for the other parameters
    for ph2 = 2:length(idx)
        tmp_struct = processed_data(idx(ph2));
        tmp_struct.runs = []; % remove the run details so the two structs can be compared for the other parameters
        % if the results are identical in terms of phase, level, design, etc then merge them
        if isequal(ref_struct,tmp_struct)
            processed_data(idx(1)).runs = [processed_data(idx(1)).runs processed_data(idx(ph2)).runs];
            phases_to_remove = [phases_to_remove idx(ph2)];
        end
    end
end
processed_data(phases_to_remove) = [];

%% compute RT per run and phase
results = struct();
for ph = 1:length(processed_data)
    results.(processed_data(ph).phase) = [];
    results.(processed_data(ph).phase).level = processed_data(ph).level;
    % initialise variables
    results.(processed_data(ph).phase).RT = [];
    results.(processed_data(ph).phase).response = [];
    
    % start processing in a run by run basis
    fprintf('=======================================\n');
    fprintf('Results for phase: "%s_%d"\n',processed_data(ph).phase,ph);
    for r = 1:length(processed_data(ph).runs)
        results.(processed_data(ph).phase).RT(r,:) = [processed_data(ph).runs(r).trials.RT];
        results.(processed_data(ph).phase).prob(r,:) = [processed_data(ph).runs(r).trials.prob];
        results.(processed_data(ph).phase).response(r,:) = [processed_data(ph).runs(r).trials.response];
        fprintf('---------------- Run %d ----------------\n',r);
        fprintf('Number of trials = %d\n', length(processed_data(ph).runs(r).trials));
        time_diff = datevec(processed_data(ph).runs(r).endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data(ph).runs(r).starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
        duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
        fprintf('Duration of run = %.2f seconds\n', duration);
        fprintf('Response Time (average) for high-probability transitions = %.2f seconds\n', nanmean(results.(processed_data(ph).phase).RT(r,results.(processed_data(ph).phase).prob(r,:)=='H' & [processed_data(ph).runs(r).trials.correct]==1))/1000);
        fprintf('Response Time (average) for low-probability transitions = %.2f seconds\n', nanmean(results.(processed_data(ph).phase).RT(r,results.(processed_data(ph).phase).prob(r,:)=='L' & [processed_data(ph).runs(r).trials.correct]==1))/1000);
        fprintf('Incorrect responses = %d\n', sum([processed_data(ph).runs(r).trials.correct]==0));
        fprintf('Missed responses = %d\n', sum(isnan([processed_data(ph).runs(r).trials.response])));
        fprintf('---------------------------------------\n');
    end
    fprintf('=======================================\n');
%     % make plots and save the figures
%     if length(processed_data(ph).runs) >= 2
%         rsp = input('Do you want to plot the Performance Index results? (y/n)','s');
%         if rsp == 'y'
%             fh1 = plot_pi(results.(processed_data(ph).phase).pi_abs,results.(processed_data(ph).phase).pi_rand,1); % 3rd argument is for fitting a curve (=1) or not (=0)
%         end
%     end
end

%% save results in the output file
if length(fname) > 1
    out_name = strcat(data.results.userref,'_',task_name);
else
    out_name = files(1:end-4);
end
save(strcat(path,out_name,'_results.mat'),'results','processed_data');
if exist('fh1','var')
    saveas(fh1,strcat(path,out_name,'_PI.fig'));
end

%% function to plot the Performance Index data per run
function figure_handle = plot_pi(pi_abs,pi_rand,fitline)
    
    % if it's not set then fit a line for PI across blocks
    if nargin < 3
        fitline = 1;
    end
    
    % plot PI points and fitted curve
    figure_handle = figure();
    hold on;
    plot(1:length(pi_abs),pi_rand,'LineWidth',1,'LineStyle','--','Color',[0.3 0.3 0.3],'DisplayName','Random guess');
    plot(1:length(pi_abs),pi_abs,'MarkerSize',30,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','.','LineWidth',1,'LineStyle','none','DisplayName','Training data','Color','k');
    if length(pi_abs) > 1 && fitline == 1
        fito = fit([1:length(pi_abs)]', pi_abs, fittype( 'a*log(x)+c', 'independent', 'x', 'dependent', 'y' ), 'StartPoint',[1 pi_abs(1)]);
        a = plot(fito);
        set(a, 'LineWidth',1,'Color','k');
    end
    
    % set axis and legend
    xlim([0 length(pi_abs)+1]);
    ylim([0 1]);
    set(gca,'FontSize',14);
    xlabel('Run','FontSize',16);
    ylabel('Performance index','FontSize',16);
    legend1 = legend('Random guess','Training data');
    set(legend1,'Location','SouthEast');
end
