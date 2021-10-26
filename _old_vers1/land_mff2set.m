function [file_name, eeg_struct] = land_mff2set(cfg, file_name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     cd(fullfile(cfg.path.project, cfg.path.data_mff))
%     folder_list = dir('*.mff')
%     %exist(folder_name, dir('*.mff'))
%     file_name = [];
%     for i_folder = 1:length(folder_list)
%         %i_folder=
%         %if contains(folder_list(i_folder).name, folder_name)
%         if contains(folder_list(i_folder).name, cfg.subj_id) 
%             if contains(folder_list(i_folder).name, cfg.session_id) 
%                 if contains(folder_list(i_folder).name, cfg.task_name,'IgnoreCase',true) 
%             
%                     file_name = folder_list(i_folder).name
%                     % no mff extension
%                     file_name = file_name(1:end-4);
%                     disp('mff file EXIST !!')
%                 end
%             end
%         end
%     end
%     if isempty(file_name); 
%         disp('!!! mff file does NOT EXIST !!'); 
%         return
%     end

    % use EEGLAB for loading the mff file 
    %cd(cfg.path.eeglab)
    %eeglab%()

    % LOAD - - - - - - 
    eeg_struct = [];
    try
        %eeg_struct = pop_mffimport({fullfile(cfg.path.project, cfg.path.data_mff, [file_name '.mff'])},{'code'});
        eeg_struct = pop_mffimport(file_name,{'code'});

    catch error_message
        disp(error_message)
        if ~isempty(error_message)
            % that means there are some problems with the triggers:
            % NO EVENTs in the data file
            %eeg_struct = pop_mffimport({fullfile(cfg.path.project, cfg.path.data_mff, [file_name '.mff'])});
            eeg_struct = pop_mffimport(file_name);
        end
    end

    % check if CHANNEL LAYOUT is correctly imported
    figure; 
    topoplot([], eeg_struct.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', ...
                 eeg_struct.chaninfo);
end

