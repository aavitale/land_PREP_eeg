function [eeg_struct, event_start_msec, event_stop_msec] = land_segment_clean(cfg, eeg_struct)
%function [eeg_struct, event_start_sampleidx, event_stop_sampleidx] = land_segment_clean(cfg, eeg_struct)
%
%   lad_segment_function 
%   takes as INPUT: eeg_struct (in eeglab format) 
%                   cfg structure: should contains the task_name (inscapes, resting_state...
%                                  and the trigger table (with the trigger name for the onset/offset events)
%
%   the function looks for the trigger name corresponding to the initial/onset event
%   / if the onset event exist ->  
%
%   OUPUT: eeg_struct "clean"

    %EEG =eeg_struct; eeglab redraw

    %!!! working only for the new folder name     
    %file_name = split(cfg.folder_name,'_')
    task_name = cfg.task_name
    
    time_vector_msec = eeg_struct.times;

    % - - - - - - - - - - - - - - - - - - - 
    % find TASK index in the trigger table
    % - - - - - - - - - - - - - - - - - - - 
    %find(strfind(task_name, trigger_table))
    %[~, task_idx] = ismember(cfg.task_name, cfg.trigger_table.task);
    % !!! if the task name does not correspond exactly to the trigger table:
    task_idx = []; 
    %trigger_id = {};
    for i_task = 1:size(cfg.trigger_table,1)
        %if contains(cfg.trigger_table.task(i_task), cfg.task_name(2:end));
        if contains(cfg.trigger_table.task(i_task), cfg.task_name, 'IgnoreCase',true);
            task_idx = [ task_idx, i_task ]
            %trigger_id{end+1} = cfg.trigger_table.trigger_start{i_task}
            trigger_start_id = cfg.trigger_table.trigger_start{i_task}
            trigger_stop_id = cfg.trigger_table.trigger_stop{i_task}
            event_duration_sec = str2num(cfg.trigger_table.duration_sec{i_task})
        
        end
    end
    if length(task_idx) > 1 
        disp('!!! problem with the trigger table')
        return
    end
    
    
    % - - - - - - - - - - - - - - - - - - - 
    % search for the TASK name in the event structure
    % - - - - - - - - - - - - - - - - - - - 
    event_start_msec = [];                     
    event_stop_msec = []; %time_vector_msec(end);  
    event_start_idx = []; event_stop_idx = [];
    
    n_event = length(eeg_struct.event);
    for i_event = 1:n_event
        if strcmp(eeg_struct.event(i_event).code, trigger_start_id)
            event_start_idx = [ event_start_idx; i_event ];
            event_start_msec = [ event_start_msec; eeg_struct(i_event).event.latency ];
  
        elseif strcmp(eeg_struct.event(i_event).code, trigger_stop_id)
            event_stop_idx = [ event_stop_idx; i_event ];
            event_stop_msec = [ event_stop_msec; eeg_struct(i_event).event.latency ];
            
        end
    end
    disp([ event_start_msec  event_stop_msec ])
    
    % check LENGTH of the acquisition
    %if ~isempty(event_stop_msec) 
    if time_vector_msec(end) < (event_start_msec + event_duration_sec*1000) 
         disp('!!! problem with the trigger table')
    end
    
    % option a): POP_RMDAT - - - - - - - - - - - - -
    % if THERE's NO START EVENT (and STOP EVENT):
    if isempty(event_start_idx) && isempty(event_stop_idx)
        %event_start_sampleidx = 1
        event_start_msec = 0 %<<<<< is it better to remove the ~first second(s) of the acquisition ???
        eeg_struct = pop_select( eeg_struct, 'time',[ event_start_msec  (event_start_msec/1000+event_duration_sec) ] );

    % if there's NO STOP trigger than cut the file according to the duration (in sec) of the task
    elseif ~isempty(event_start_idx) && isempty(event_stop_idx)
        eeg_struct = pop_rmdat( eeg_struct, {trigger_start_id}, ...
                        [ 0  event_duration_sec ] ,0);

    elseif ~isempty(event_start_idx) && ~isempty(event_stop_idx)
        eeg_struct = pop_rmdat( eeg_struct, {trigger_start_id}, ...
                        [ event_start_msec/1000   event_stop_msec/1000 ] ,0);
        
    % else isempty(event_start_idx) && ~isempty(event_stop_idx)  % ?????
    
    end
    
    
%     % option b): POP_SELECT !!! less precise - - - - - - - - - - - - -
%     % identify the sample corresponding to the trigger_start and remove the samples before it
%     % it could generate some errors
% 
%     %event_name = {}; event_latency = [];
%     event_start_msec = []; event_stop_msec = [];
%     event_start_sampleidx = []; event_stop_sampleidx = [];
%     for i_event = 1:n_event
%         if strcmp(eeg_struct(i_event).event.code, cfg.trigger_table.trigger_start{task_idx})
%         %if contains(cfg.trigger_table.trigger_start{task_idx}, eeg_struct(i_event).event.code)
%             event_start_msec = eeg_struct(i_event).event.latency
%             %event_name{i_event,1} = eeg_struct(i_event).event.code;
%             %event_latency[i_event,1] = eeg_struct(i_event).event.latency;
% 
%             %[ ~, event_start_sampleidx ] = min(abs(time_vector_msec - (event_start_msec)));
%             [ ~, event_start_sampleidx ] = min(abs(time_vector_msec - floor(event_start_msec)));
% 
%         elseif strcmp(eeg_struct(i_event).event.code, cfg.trigger_table.trigger_stop{task_idx})
%         %elseif contains(eeg_struct(i_event).event.code, cfg.trigger_table.trigger_stop{task_idx})
%             event_stop_msec = eeg_struct(i_event).event.latency
%             [ ~, event_stop_sampleidx ] = min(abs(time_vector_msec - event_stop_msec))
% 
%         end    
%     end
% 
%     if isempty(event_start_sampleidx)
%         event_start_sampleidx = 1
%     end
%     
%     % if there's no STOP trigger than cut the file according to the duration (in sec) of the task
%     if isempty(event_stop_sampleidx)
%         %event_stop_sampleidx = length(time_vector_msec) -1
%         
%         video_duration_msec = str2num(cfg.trigger_table.duration_sec{task_idx})*1000;
%         %[ ~, event_stop_sampleidx ] = min(abs(time_vector_msec - video_duration_sec*1000));
%         [ ~, event_stop_sampleidx ] = min(abs(time_vector_msec - ...
%                                                 (event_start_sampleidx + video_duration_msec -1)));
%     end
% 
%     %eeg_struct = pop_select( eeg_struct, 'point',[event_start_sampleidx  event_stop_sampleidx+1] );
%     eeg_struct = pop_select( eeg_struct, 'point',[event_start_sampleidx-1  event_stop_sampleidx+1] );
    
    % - - - - - - - - -  -
    % SCROLL the data to check that the initial trigger is at sample = 1
    pop_eegplot(eeg_struct, 1, 1, 1);
    % EEG = eeg_struct; eeglab redraw
    
    % BOUNDARY and EVENT should belong to 2 different samples
    if length(eeg_struct.event) == 0
        disp('!!! no event trigger in the file')
    else
        eeg_struct.event(1).latency
        eeg_struct.event(2).latency
        if floor(eeg_struct.event(2).latency) == floor(eeg_struct.event(1).latency)
            disp('!!! problem with the boundary position -> downsampling is going to crash')
            
            eeg_struct.event(2).latency = event_start_msec/1000
            %eeg_struct.event(1).duration = []; )
        
        end
    end
end
