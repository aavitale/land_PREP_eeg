function [file_name] = land_findfile(cfg, folder_name, file_ext)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %cd(fullfile(cfg.path.project, cfg.path.data_mff))
    %folder_list = dir('*.mff')
    folder_list = dir(['*.' file_ext])
    %exist(folder_name, dir('*.mff'))
    file_name = [];
    for i_folder = 1:length(folder_list)
        %i_folder=
        %if contains(folder_list(i_folder).name, folder_name)
        if contains(folder_list(i_folder).name, cfg.subj_id) 
            if contains(folder_list(i_folder).name, cfg.session_id) 
                if contains(folder_list(i_folder).name, cfg.task_name,'IgnoreCase',true) 
            
                    file_name = folder_list(i_folder).name
                    % no mff extension ?
                    %file_name = file_name(1:end-4);
                    disp([ file_ext ' file EXIST !!'])
                end
            end
        end
    end
    if isempty(file_name); 
        disp('!!! file does NOT EXIST !!'); 
        return
    end
end
