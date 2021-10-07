%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 25/03/2020
% This script analyses the Questionnaire task results collected from the iABC app
% In order for this to work, it requires that the read_iABC_results script has been run prior to this
%% clear workspace and figures
clearvars; clc; close all;
%% read results file
task_name = 'qn';
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
%% load the file and printout the answers
load(fname{1});
for i = 1:length(results.answers)
    if iscell(results.answers{i}.answer)
        fprintf('Answer to Q%d: %s\n',i,results.answers{i}.answer{1});
    else
        fprintf('Answer to Q%d: %s\n',i,results.answers{i}.answer);
    end
end
