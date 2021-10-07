%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 30/01/2020
% This script analyses the Sequence Learning task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this

%% clear workspace and figures
% clearvars; clc; close all;
%% read results file
task_name = 'SL';
ignore_practice = 1; % This will skip the results for the practice phase (=1). Turn this to =0 if you want them in the output.
remove_short_runs = 30; % Set this number to the minimum acceptable trials (shorter runs will be removed!)
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
            processed_data(phase_idx).symbolset = data.results.items{tr}.symbolset;
        end
    end
    % find the level and design information
    if isfield(data.results.items{tr},'condition') % revisit this if multiple conditions are allowed per phase!!!
        processed_data(phase_idx).level = data.results.items{tr}.level;
        processed_data(phase_idx).symbolMapping = data.results.items{tr}.symbolMapping;
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
            processed_data(phase_idx).runs(block_idx).blockSequence = data.results.items{tr}.blockSequence' + 1;
            processed_data(phase_idx).runs(block_idx).starttime = data.results.items{tr}.start;
        else
            processed_data(phase_idx).runs(block_idx).endtime = data.results.items{tr}.end;
            processed_data(phase_idx).runs(block_idx).PI = data.results.items{tr}.blockPerformanceIndex;
            processed_data(phase_idx).runs(block_idx).seq_dist = [];
            processed_data(phase_idx).runs(block_idx).resp_dist = [];
            switch processed_data(phase_idx).level
                % for level-0 just read the saved distributions
                case 0
                    processed_data(phase_idx).runs(block_idx).seq_dist = str2num(data.results.items{tr}.sequenceDistribution);
                    processed_data(phase_idx).runs(block_idx).resp_dist = str2num(data.results.items{tr}.responseDistribution);
                % for level-1 convert it into a 4x4 matrix
                case 1
                    tmp = strsplit(data.results.items{tr}.sequenceDistribution,{'[',']'});
                    for i = 1:length(tmp)
                        if ~isempty(tmp{i})
                            processed_data(phase_idx).runs(block_idx).seq_dist = [processed_data(phase_idx).runs(block_idx).seq_dist; str2num(tmp{i})];
                        end
                    end
                    tmp = strsplit(data.results.items{tr}.responseDistribution,{'[',']'});
                    for i = 1:length(tmp)
                        if ~isempty(tmp{i})
                            processed_data(phase_idx).runs(block_idx).resp_dist = [processed_data(phase_idx).runs(block_idx).resp_dist; str2num(tmp{i})];
                        end
                    end
                % for level-2 convert it into a 16x4 matrix
                case 2
                    tmp = strsplit(data.results.items{tr}.sequenceDistribution,{'[',']'});
                    for i = 1:length(tmp)
                        if ~isempty(tmp{i})
                            processed_data(phase_idx).runs(block_idx).seq_dist = [processed_data(phase_idx).runs(block_idx).seq_dist; str2num(tmp{i})];
                        end
                    end
                    tmp = strsplit(data.results.items{tr}.responseDistribution,{'[',']'});
                    for i = 1:length(tmp)
                        if ~isempty(tmp{i})
                            processed_data(phase_idx).runs(block_idx).resp_dist = [processed_data(phase_idx).runs(block_idx).resp_dist; str2num(tmp{i})];
                        end
                    end
            end
        end
    end
    % find the information per trial (if PI and distributions above are correct, then you might not need this!!!)
    if isfield(data.results.items{tr},'trial')
        trial_idx = trial_idx + 1;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).sequenceItems = data.results.items{tr}.sequence + 1;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).gapSize = data.results.items{tr}.gapSize;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).targetPos = data.results.items{tr}.targetPos;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).drawOrder = str2num(data.results.items{tr}.drawOrder) + 1;
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).correctAnswer = data.results.items{tr}.correctAnswer + 1;
        if isfield(data.results.items{tr},'response')
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).response = data.results.items{tr}.response + 1;
        else
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).response = NaN;
        end
        processed_data(phase_idx).runs(block_idx).trials(trial_idx).RT = data.results.items{tr}.responseTime;
        if sum(data.results.items{tr}.durationJitter~=0) > 0
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).tJitter = data.results.items{tr}.durationJitter;
        end
        if sum(data.results.items{tr}.jitterX~=0) > 0
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).sJitterX = data.results.items{tr}.jitterX;
        end
        if sum(data.results.items{tr}.jitterY~=0) > 0
            processed_data(phase_idx).runs(block_idx).trials(trial_idx).sJitterY = data.results.items{tr}.jitterY;
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

%% compute performance and strategy per run and phase
results = struct();
for ph = 1:length(processed_data)
    results.(processed_data(ph).phase) = [];
    results.(processed_data(ph).phase).level = processed_data(ph).level;
    % initialise variables
    results.(processed_data(ph).phase).RT = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).pi_abs = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).pi_rand = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).pi_rel = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).strategy = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).exactMatch = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).exactMaximize = nan(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).seq_distr = cell(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).resp_distr = cell(length(processed_data(ph).runs),1);
    results.(processed_data(ph).phase).maxi_distr = cell(length(processed_data(ph).runs),1);
    % start processing in a run by run basis
    fprintf('=======================================\n');
    fprintf('Results for phase: "%s_%d"\n',processed_data(ph).phase,ph);
    for r = 1:length(processed_data(ph).runs)
        fprintf('---------------- Run %d ----------------\n',r);
        fprintf('Number of trials = %d\n', length(processed_data(ph).runs(r).trials));
        time_diff = datevec(processed_data(ph).runs(r).endtime,'yyyy-mm-ddTHH:MM:SS') - datevec(processed_data(ph).runs(r).starttime,'yyyy-mm-ddTHH:MM:SS');    % decimal points for seconds are removed due to precision limit
        duration = time_diff(1)*365*24*60*60 + time_diff(2)*30*24*60*60 + time_diff(3)*24*60*60 + time_diff(4)*60*60 + time_diff(5)*60 + time_diff(6);	% 4th column is hours, 5th is minutes, 6th is seconds
        fprintf('Duration of run = %.2f seconds\n', duration);
        results.(processed_data(ph).phase).RT(r) = mean([processed_data(ph).runs(r).trials(isfinite([processed_data(ph).runs(r).trials.response]')).RT])/1000;
        fprintf('Response Time (average) = %.2f seconds\n', results.(processed_data(ph).phase).RT(r));
        fprintf('Missed responses = %d\n', sum(isnan([processed_data(ph).runs(r).trials.response])));
        % calculate PI
        [results.(processed_data(ph).phase).pi_abs(r), results.(processed_data(ph).phase).pi_rand(r), results.(processed_data(ph).phase).pi_rel(r), results.(processed_data(ph).phase).seq_distr{r}, results.(processed_data(ph).phase).resp_distr{r}] = calculate_pi([processed_data(ph).runs(r).trials.correctAnswer]', [processed_data(ph).runs(r).trials.response]', {processed_data(ph).runs(r).trials.sequenceItems}', processed_data(ph).level, size(processed_data(ph).design));
        fprintf('Performance Index (absolute) = %2.2f%%\n', 100*results.(processed_data(ph).phase).pi_abs(r));
        % calculate Strategy choice
        [results.(processed_data(ph).phase).strategy(r), results.(processed_data(ph).phase).exactMatch(r), results.(processed_data(ph).phase).exactMaximize(r), results.(processed_data(ph).phase).maxi_distr{r}] = calculate_strategy([processed_data(ph).runs(r).trials.correctAnswer]', [processed_data(ph).runs(r).trials.response]', {processed_data(ph).runs(r).trials.sequenceItems}', processed_data(ph).level, processed_data(ph).design);
        fprintf('Strategy choice = %1.3f\n', results.(processed_data(ph).phase).strategy(r));
        fprintf('---------------------------------------\n');
    end
    % calculate Strategy ICD if runs are more than 1
    if length(processed_data(ph).runs) >= 2
        [results.(processed_data(ph).phase).strategy_icd, results.(processed_data(ph).phase).strategy_ini, results.(processed_data(ph).phase).strategy_end, results.(processed_data(ph).phase).strategy_change] = calculate_strategy_icd(results.(processed_data(ph).phase).strategy,results.(processed_data(ph).phase).exactMatch);
        fprintf('Strategy ICD for phase "%s_%d" = %1.3f\n', processed_data(ph).phase,ph, results.(processed_data(ph).phase).strategy_icd);
    else
        fprintf('Strategy ICD can''t be defined for this phase.\n');
    end
    fprintf('=======================================\n');
    % make plots and save the figures
    if length(processed_data(ph).runs) >= 2
        rsp = input('Do you want to plot the Performance Index results? (y/n)','s');
        if rsp == 'y'
            fh1 = plot_pi(results.(processed_data(ph).phase).pi_abs,results.(processed_data(ph).phase).pi_rand,1); % 3rd argument is for fitting a curve (=1) or not (=0)
        end
        rsp = input('Do you want to plot the Strategy results? (y/n)','s');
        if rsp == 'y'
            fh2 = plot_strategy(results.(processed_data(ph).phase).strategy,results.(processed_data(ph).phase).exactMatch,results.(processed_data(ph).phase).exactMaximize);
        end
    end
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
if exist('fh2','var')
    saveas(fh2,strcat(path,out_name,'_Strategy.fig'));
end

%% function to calculate Performance Index from the input data (it should be called per run)
function [pi_abs, pi_rand, pi_rel, seq_distr, resp_distr] = calculate_pi(answers, response, items, level, dim)
    
    % get context per trial
    if level > 0
        context = zeros(length(answers),level);
    end
    for tr = 1:length(answers)
        context(tr,:) = items{tr}(end-level+1:end);
    end
    
    % calculate response and sequence distributions
    resp_distr = zeros(dim);
    seq_distr = zeros(dim);
    
    for tr = 1:size(context,1)
        
        switch level
        case 0              % for level-0 and random sequence, we have no context
            context_idx = 1;
        case 1              % for level-1 context is the last symbol of the sequence
            context_idx = context(tr,end);
        case 2              % for level-2 context is the last two symbols but in the following order (e.g. {11},{21},{31},{41},{12},{22},{32},{42},etc)
            context_idx = context(tr,end-1) + (context(tr,end) - 1) * dim(2);
        end
        
        seq_distr(context_idx,answers(tr)) = seq_distr(context_idx,answers(tr)) + 1;   % increase counter

        if isnan(response(tr))	% for no response
            resp_distr(context_idx,:) =  resp_distr(context_idx,:) + 0.25;
        else
            resp_distr(context_idx,response(tr)) = resp_distr(context_idx,response(tr)) + 1;
        end
    end
    
    % calculate PI
    pi_abs = sum( sum( min(seq_distr,resp_distr),2) .* 1 );    % only one PI per session
    pi_rand = sum( sum( min((sum(seq_distr,2)/size(seq_distr,2))*ones(1,size(seq_distr,2)),seq_distr),2) .* 1 );    % sum along targets for random
    
    pi_abs = pi_abs / length(answers);      % normalise by number of trials
    pi_rand = pi_rand / length(answers);    % normalise by number of trials
    pi_rel = pi_abs - pi_rand;              % PI relative (= absolute - random)
    
end

%% function to calculate Strategy ICD from the input data (it should be called per run)
function [strategy, exactMatch, exactMaximize, maxi_distr] = calculate_strategy(answers, response, items, level, model_distr)
    
    % get context per trial
    if level > 0
        context = zeros(length(answers),level);
    end
    for tr = 1:length(answers)
        context(tr,:) = items{tr}(end-level+1:end);
    end
    
    % calculate response and sequence distributions
    resp_distr = zeros(size(model_distr));
    seq_distr = zeros(size(model_distr));
    
    for tr = 1:size(context,1)
        
        switch level
        case 0              % for level-0 and random sequence, we have no context
            context_idx = 1;
        case 1              % for level-1 context is the last symbol of the sequence
            context_idx = context(tr,end);
        case 2              % for level-2 context is the last two symbols but in the following order (e.g. {11},{21},{31},{41},{12},{22},{32},{42},etc)
            context_idx = context(tr,end-1) + (context(tr,end) - 1) * size(model_distr,2);
        end
        
        seq_distr(context_idx,answers(tr)) = seq_distr(context_idx,answers(tr)) + 1;   % increase counter
        
        if isnan(response(tr))	% for no response
            resp_distr(context_idx,:) =  resp_distr(context_idx,:) + 0.25;
        else
            resp_distr(context_idx,response(tr)) = resp_distr(context_idx,response(tr)) + 1;
        end
    end
    
    % create maximization distribution
    ideal = model_distr * length(answers);                          % get ideal model counts for each context-target combination
    maxi_distr = zeros(size(model_distr));                          % initialise maximization matrix
    
    for i = 1:size(maxi_distr,1)                                    % for each context, do
        maxItem = find(ideal(i,:) == max(ideal(i,:)));              % find the target with the maximum count
        maxi_distr(i, maxItem) = sum(ideal(i,:));                   % update maximization matrix for the above target and give it the sum of all targets' counts
    end
    
    % renormalise distribution matrices
    seq_matrix = seq_distr + 0.1;                       % we add 0.1 because we are going to use Kullback-Leibler divergence which requires non zero values
    resp_matrix = resp_distr + 0.1;
    maxi_matrix = maxi_distr + 0.1;
    
    seq_matrix = seq_matrix/sum(sum(seq_matrix));       % normalise after adding 0.1
    resp_matrix = resp_matrix/sum(sum(resp_matrix));
    maxi_matrix = maxi_matrix/sum(sum(maxi_matrix));
    
    seq_matrix = seq_matrix(:);                         % converts an (NxM) matrix to an (NxM) vector for easier computations
    resp_matrix = resp_matrix(:);
    maxi_matrix = maxi_matrix(:);
    
    % calculate strategy
    Relative = 0;
    Maximize = 0;
    ExtMatch = 0;
    ExtMaximize = 0;
    for i = 1:length(seq_matrix)
        Relative = Relative + seq_matrix(i) * log(seq_matrix(i)/resp_matrix(i));       % KL divergence
        Maximize = Maximize + maxi_matrix(i) * log(maxi_matrix(i)/resp_matrix(i));
        ExtMatch = ExtMatch + maxi_matrix(i) * log(maxi_matrix(i)/seq_matrix(i));
        ExtMaximize = ExtMaximize + seq_matrix(i) * log(seq_matrix(i)/maxi_matrix(i));
    end
    strategy = Maximize - Relative;
    exactMatch = ExtMatch;
    exactMaximize = -ExtMaximize;
    
end

%% function to calculate Strategy ICD from strategy and exactMatch values
function [strategy_icd, strategy_ini, strategy_end, strategy_change] = calculate_strategy_icd(all_strategy,all_exactMatch)

    area_match = trapz(1:length(all_exactMatch), all_exactMatch);
    area_resp  = trapz(1:length(all_strategy), all_strategy);
    
    strategy_icd = (area_match - area_resp)/length(all_strategy);                       % difference between response and matching curves
    strategy_end = mean(all_exactMatch(end-1:end)) - mean(all_strategy(end-1:end));     % end strategy is the difference of the last two runs
    strategy_ini = mean(all_exactMatch(1:2)) - mean(all_strategy(1:2));                 % ini strategy is the difference of the first two runs
    strategy_change = strategy_end - strategy_ini; 

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

%% function to plot the Strategy data per run
function figure_handle = plot_strategy(strategy,exactMatch,exactMaximize)
    
    % plot Strategy curves
    figure_handle = figure();
    hold on;
    plot(1:length(strategy), strategy, 'MarkerSize',10,'Marker','o','LineWidth',2,'Color',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3]);
    plot(1:length(strategy), exactMatch,'LineStyle','--','LineWidth',3,'Color',[0 0 0],'Marker','none');
    plot(1:length(strategy), exactMaximize,'LineWidth',3,'Color',[0 0 0],'Marker','none');
    
    % set axis and legend
    xlim([0 length(strategy)+1]);
    ylim([-1 1]);
    set(gca,'FontSize',14);
    xlabel('Run','FontSize',16);
    ylabel('{\Delta}KL(Model_M_a_x_i_m_i_s_a_t_i_o_n,Model_M_a_t_c_h_i_n_g)','FontSize',16);
    legend1 = legend('Strategy','Exact Matching','Exact Maximisation');
    set(legend1,'Location','NorthEast');
    
end
