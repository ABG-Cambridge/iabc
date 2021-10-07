%% Script written by Vasilis Karlaftis (vmk25@cam.ac.uk) 30/01/2020
% This script reads output data from the iABC app in json format and
% converts them into a Matlab struct for further analysis

%% select and read files using GUI (multiple files are allowed to be selected)
[files,path,msg] = uigetfile('*.json','MultiSelect','on');
if msg == 0
    fprintf('* No file was selected *\n');
    return;
end
fname = strcat(path,files);
if ischar(fname) == 1
    fname = {fname};
end
%% repeat the following for each selected file
for f = 1:length(fname)
    % open file for reading
    fid = fopen(fname{f});
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    % convert file data into a struct
    data = jsondecode(str);
    
    %% for task in the output file (if all data are downloaded for a single subject)
    for test = 1:length(data)
        % select data appropriately based on their format
        if isstruct(data)
            testdata = data(test);
        else
            testdata = data{test};
        end
        % for questionnaire task, store task details until reached the answers (stored at the end of the file)
        if isfield(testdata,'taskId')
            if strcmp(testdata.taskId,'qn')
                store_data = testdata;
                continue;
            end
        else
            % for questionnaire answers, append them on the task info (it only works if they are saved in this order)
            fields = fieldnames(store_data);
            for fd = 1:length(fields)
                testdata.(fields{fd}) = store_data.(fields{fd});
            end
        end
        % make a printout of the file info (e.g. task, app version, etc)
        fprintf('----- Results are being processed -----\n');
        fprintf('User: %s\n',testdata.userref);
        fprintf('App Version: %s\n',testdata.sessionInfo.appVersion);
%         if isstruct(testdata.items)
%             fprintf('Project: %s\n',testdata.items(1,1).project);
%         else
%             fprintf('Project: %s\n',testdata.items{1,1}.project);
%         end
        fprintf('Task: %s\n',testdata.taskId);
        fprintf('Task state: %s\n',testdata.state);
        % if state is not "completed", ask user if wishes to save the output
        if ~strcmp(testdata.state,'completed')
            txt = strcat({'Task state is '},testdata.state,{'. Do you want to save these data? (y/n)'});
            rsp = input(txt{1},'s');
            if rsp == 'n'
                continue;
            end
        end
        
        % save data into a matlab file for further analysis
        outname = strcat(path,testdata.userref,'_',testdata.taskId,'.mat');
        ext = 1;
        while exist(outname,'file')
            outname = strcat(path,testdata.userref,'_',testdata.taskId,'_',num2str(ext),'.mat');
            ext = ext + 1;
        end
        results = testdata;
        save(outname,'results');
        fprintf('---------------- Done -----------------\n');
        
    end
    
end
