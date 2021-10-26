function eeg_prep = land_PREprocess_run(cfg, eeg_struct) 
% 
% preprocess LAND (segmented) RestingState or VIdeos data
% with LAND preprocessing pipeline
%
%
%   INPUT
%       datafile = full filename to raw data *.set file
%       !!! TO IMPLEMENT: run_on_server = true or false
%
%   Example usage:
%
%
%   SAVE data  at 4 different stages:
%   () eye-struct
%   - after band-pass + NOTCH filtering
%   - after ASR 
%           vis_artifacts(eeg_prep.eeg_afterASR, eeg_prep.eeg_notch)
%   - after ICA decomposition (without badICs removal)
%   - after badICs removal
%           
%  VERSION: - - - - - - - - - -
%       20211022: channel selection (inteprolation) based on ASR + HURST 
%                   to reduce the number of channels interpolated
%
%       20211025: remove the channel selection (inteprolation) 
%                   based on ASR + HURST -> leave only ASR%       


%%
% function parameter_dump4report(cfg)
% 
% fprintf('Subject ID: %s \n',cfg.subid);
% fprintf('Task: %s \n',cfg.dataType);
% fprintf('File: %s \n',cfg.setfilename);
% fprintf('\n');
% fprintf('Downsampling rate: %d \n',cfg.pp.downsample_rate);
% fprintf('High pass cutoff: %d \n',cfg.pp.hpf_cutoff);
% fprintf('Low pass cutoff: %d \n',cfg.pp.lpf_cutoff);
% fprintf('Line noise notch frequency: %d Hz + or - %d Hz \n', ...
%     cfg.pp.line_noise_freq, ...
%     cfg.pp.line_noise_boundaries);
% fprintf('Outer channels to remove: \n');
% cfg.pp.chan_toreject
% 
% if cfg.pp.do_chan_interp_pruning
%     fprintf('Interpolate bad channels: true \n')
% else
%     fprintf('Interpolate bad channels: false \n')
% end % if cfg.pp.do_chan_interp_pruning
% 
% fprintf('Channels to use for interpolation and then prune before ICA: \n');
% cfg.pp.chan_interp_prune
% 
% fprintf('IClabels parameters to use when classifying components \n');
% cfg.pp.iclabels_params
% 
% fprintf('ASR parameters: \n')
% cfg.pp.asr
% 
% end % parameter_dump4report


% DATASET already imported/converted and segmented (excluding bad portions)
% fprintf('Loading dataset \n')
% eeg_struct = pop_loadset('filename',cfg.setfilename);
% eeg_raw = eeg_struct;


%% Checks
sample_rate = eeg_struct.srate;
n_sample = eeg_struct.pnts;
n_chan = eeg_struct.nbchan;
n_chan_max = sqrt(n_sample/20);
% if n_sample > n_chan^2 * 20
%     disp([ num2str ' channels can be given as input to ICA'])
% else
%     sprintf('Number of channels for ICA should be reduced to %d' n_chan_max);
% end


% ####################################################################
%% Subset of 11 CHANNELS as EOG
% ####################################################################
eye_chan = cfg.prep.eye_channels;
eye_struct = pop_select(eeg_struct, 'channel', eye_chan);
%%pop_chanedit(eye_struct.chanlocs)

% if cfg.do_save_eye
%     eye_struct = pop_saveset( eye_struct, 'filename', [file_name '_EOG.set'], ...
%          'filepath', cfg.path.data_set)
%      % or in the eeg_prep general structure ??
% end


% ####################################################################
%% REMOVE external ring of channels
% ####################################################################
fprintf('Channel removal \n')
chan_toreject = cfg.prep.chan_toreject;
eeg_struct = pop_select(eeg_struct, 'nochannel', chan_toreject);
eeg_raw_chanred = eeg_struct;
% eeg_struct = eeg_raw_chanred;


%% ####################################################################
% FILTERING
% ####################################################################
%% Downsampling
% fprintf('Downsampling \n')
% eeg_struct = pop_resample(eeg_struct, cfg.prep.downsample_rate);
% % fix latency so they are integers
% % for i = 1:length(eeg_struct.event)
% %     eeg_struct.event(i).latency = ceil(eeg_struct.event(i).latency);
% % end % for i
% eeg_down = eeg_struct;
% 
% 
%% HIGH-pass filtering
% fprintf('Band-pass filtering \n')
% eeg_struct = pop_eegfiltnew(eeg_struct, cfg.prep.hpf_cutoff, [], [],0,[],0);
% eeg_hpf = eeg_struct;


%% LOW-PASS Filtering
eeg_struct = pop_eegfiltnew(eeg_struct, [], cfg.prep.lpf_cutoff, [],0,[],0);
eeg_lpf = eeg_struct;


%% Remove line noise
fprintf('Remove line noise with notch filter \n')
eeg_notch = pop_eegfiltnew(eeg_struct, ...
    'locutoff',cfg.prep.line_noise_freq-cfg.prep.line_noise_boundaries, ...
    'hicutoff',cfg.prep.line_noise_freq+cfg.prep.line_noise_boundaries, ...
    'revfilt',1,'plotfreqz',1);

pop_eegplot(eeg_notch, 1, 1, 1);
% if cfg.do_save_fig
%     % fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
%     %     sprintf('%s_%s_preproc_badSamplesChannels.jpg',subid,dataType));
%     fname2save = fullfile(cfg.path.data_prep, ...
%         sprintf('%s_preproc_notchdata_scroll',cfg.folder_name));
%     %print(gcf,fname2save, '-djpeg','-noui');
%     print(gcf,fname2save, '-dtiffn','-noui');    
% end

%EEG = eeg_notch; eeglab redraw

% ####################################################################
%% Iterative of BAD CHANnel and sample detection 
% ####################################################################
% with cleanraw and ASR - - - - - - - - - - 
% 1° step: DETECT BAD CHANNELS (and samples) with ASR
% 2° step: detect BAD CHAN with FASTER metrics
% 3° step: BAD CHAN intersect (NO!!!)  

% 4° step: detect remaining BAD sample (with ASR)

% 5° step: REPAIR the timeseries without removing any channel or sample

% 1° - - - - - - - - - - - - - - - - - -  - -- - -
n_iter = cfg.prep.asr.iterations;
fprintf('Iterative of bad channel and sample detection \n')

% set burst rejection = ON -> to identify also bad samples (rejected and not repaired
bad_chan_table = []; bad_sample_table = [];
for i_iter = 1:n_iter
    sprintf('Iteration %d',i_iter);

    [eeg_cleanraw_iter, HP, BUR, bad_chan_cleanraw] = clean_artifacts(eeg_notch, ...
        'ChannelCriterion',cfg.prep.asr.chancrit, ...
        'LineNoiseCriterion',cfg.prep.asr.lncrit, ...
        'FlatlineCriterion',cfg.prep.asr.flcrit, ...
        'BurstCriterion',cfg.prep.asr.burstcrit, ...
        'WindowCriterion',cfg.prep.asr.windowcrit, ...
        'Highpass',cfg.prep.asr.highpass, ...
        'BurstRejection','on');  %<<<<<<<<<<<<<<<<<<<

    bad_chan_table(i_iter,:) = ~eeg_cleanraw_iter.etc.clean_channel_mask;
    bad_sample_table(i_iter,:) = ~eeg_cleanraw_iter.etc.clean_sample_mask;
end % for i_iter = 1:n_iter

bad_sample_ASR = find(sum(bad_sample_table,1) >= cfg.prep.asr.niter_bad_chan);  % bad sample identified in more than n iterations
bad_sample_ASRperc = sum(bad_sample_table,2) ./ eeg_notch.pnts *100;

bad_chan_ASR = find(sum(bad_chan_table,1) >= cfg.prep.asr.niter_bad_chan);  % bad channel identified in more than n iterations

vis_artifacts(eeg_cleanraw_iter, eeg_notch)
%print(gcf,fname2save, '-djpeg','-noui');
if cfg.do_save_fig
    fname2save = fullfile(cfg.path.data_prep, ...
        sprintf('%s_ASRbadsample_reject',cfg.folder_name));
    %print(gcf,fname2save, '-djpeg','-noui');
    print(gcf,fname2save, '-dtiffn','-noui');    
end

% 2° -- - - - - - - - - - - - -  -
% identify BAD CHANNELS based on mean_correlation, variance and HURST exponent (as in FASTER
% as in FASTER toolbox (with a z_threshold = 2.5 
[ bad_chan_zthresh3, chan_properties, z_score ] = land_BADchan_detection_continuous(eeg_notch, 3);
[ bad_chan_zthresh2, chan_properties, z_score ] = land_BADchan_detection_continuous(eeg_notch, 2.5);
length(bad_chan_zthresh2)

bad_chan_FAST = bad_chan_zthresh2;

bad_chan_FASTidx = zeros([1,size(chan_properties,1)]);
bad_chan_FASTidx(bad_chan_zthresh2) = +10;
bad_chan_FASTidx(bad_chan_zthresh3) = +10;

% 3° -- - - - - - - - - - - - -  -
if cfg.prep.asr_fast_intersect  %<< detect bad channel based on both method 
%     % COMPARE and MERGE the results of the 2 methods:
%     %bad_chan_idx = intersect(find(bad_chan_ASR ==1), find(bad_chan_FAST ==1));
%     bad_chan_idx = intersect(bad_chan_ASR, bad_chan_FAST)
% 
%     bad_chan_label = {};
%     counter = 1;
%     for i_chan = 1:length(bad_chan_idx)
%         bad_chan_label{1,counter} = eeg_notch.chanlocs(bad_chan_idx(i_chan)).labels;
%         counter = counter+1;
%     end % for i_chan = 1:length(eeg_cleanraw.etc.clean_channel_mask)
% 
%     eeg_bad_chan = pop_select( eeg_notch, 'channel',bad_chan_label);
% 
%     % 3bis: REMOVE BAD chan:
%     eeg_nobad_chan = pop_select( eeg_notch, 'nochannel', bad_chan_label);
% 
%     %4° - - - - - - - -  -- - - -  -
%     eeg_cleanraw  = clean_artifacts(eeg_nobad_chan, ...
%             'FlatlineCriterion','off', ...%cfg.prep.asr.flcrit, ...        
%             'ChannelCriterion','off', ... %cfg.prep.asr.chancrit, ...
%             'LineNoiseCriterion','off', ... %cfg.prep.asr.lncrit, ...
%             'BurstCriterion',cfg.prep.asr.burstcrit, ...
%             'WindowCriterion','off', ...
%             'Highpass','off', ...
%             'BurstRejection','off') %<<<<<<<<<<<<<<<<<<<
% 
%     %vis_artifacts(eeg_cleanraw, eeg_nobad_chan)

    
% = = = = = = = = = = = = = = = = = = = = = =
elseif cfg.prep.asr_fast_intersect == 0
% = = = = = = = = = = = = = = = = = = = = = =
    bad_chan_idx = bad_chan_ASR;

    bad_chan_label = {};
    counter = 1;
    for i_chan = 1:length(bad_chan_idx)
        bad_chan_label{1,counter} = eeg_notch.chanlocs(bad_chan_idx(i_chan)).labels;
        counter = counter+1;
    end % for i_chan = 1:length(eeg_cleanraw.etc.clean_channel_mask)

    eeg_bad_chan = pop_select( eeg_notch, 'channel',bad_chan_label);
    
    % final attempt to REPAIR the bad samples %without removing any part of the signal
    [eeg_cleanraw, HP, BUR, bad_chan_cleanraw] = clean_artifacts(eeg_notch, ...
        'ChannelCriterion',cfg.prep.asr.chancrit, ...
        'LineNoiseCriterion',cfg.prep.asr.lncrit, ...
        'FlatlineCriterion',cfg.prep.asr.flcrit, ...
        'BurstCriterion',cfg.prep.asr.burstcrit, ...
        'WindowCriterion',cfg.prep.asr.windowcrit, ...
        'Highpass',cfg.prep.asr.highpass, ...
        'BurstRejection',cfg.prep.asr.burstrejection);  
 
    vis_artifacts(eeg_cleanraw, eeg_notch)
    %vis_artifacts(eeg_prep.eeg_afterASR, eeg_prep.eeg_notch)
    
    if cfg.do_save_fig
        fname2save = fullfile(cfg.path.data_prep, ...
            sprintf('%s_ASRbadsample_repair',cfg.folder_name));
        print(gcf,fname2save, '-dtiffn','-noui');    
    end
end

    % - - - - - - - -  - - - - - - - - - - -
    % check which portion of the data is still considered as artifactual (even after ASR)
    [eeg_cleanraw_ASR2, HP, BUR, bad_chan_cleanraw] = clean_artifacts(eeg_cleanraw, ...
        'ChannelCriterion',cfg.prep.asr.chancrit, ...
        'LineNoiseCriterion',cfg.prep.asr.lncrit, ...
        'FlatlineCriterion',cfg.prep.asr.flcrit, ...
        'BurstCriterion',cfg.prep.asr.burstcrit, ...
        'WindowCriterion',cfg.prep.asr.windowcrit, ...
        'Highpass',cfg.prep.asr.highpass, ...
        'BurstRejection','on');

% STORE this information into the cfg.prep structure:
cfg.prep.bad_chan_ASR = bad_chan_ASR;
cfg.prep.bad_chan_FAST = bad_chan_FAST;
cfg.prep.bad_chan_idx = bad_chan_idx;
cfg.prep.bad_sample_ASR = bad_sample_ASR;
cfg.prep.bad_sample_ASRperc = bad_sample_ASRperc
cfg.prep.bad_sample_ASR2 = ~eeg_cleanraw_ASR2.etc.clean_sample_mask;
    
% = = = = = = = = = = == = = = = = = = =
% TO DO ??? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% neighbouring structure
    
    
%% FIGURE plotting bad channels and samples
figure; set(gcf,'position',[10,10,1800,1000])
subplot(331) %(6,3,[1,4])
imagesc(bad_chan_table')
xlabel('iteration');
ylabel('CHANNEL')
title([ 'Bad Channels ASR = ' num2str(round(length(bad_chan_ASR)/size(bad_chan_table,2)*100)) '%)'], ...
            'FontSize', 20)

subplot(3,3,[4,7]);
topoplot(sum(bad_chan_table,1), ...
    eeg_notch.chanlocs, ...
    'electrodes','labelpoint','chaninfo', ...
    eeg_notch.chaninfo) ;
title('Bad channel ASR')
x_lim = get(gca, 'Xlim')
%text(x_lim(1)*1.5,0, bad_chan_label, 'fontweight','bold') %,'FontSize',14)

subplot(3,3,[5,8]);
topoplot(bad_chan_FASTidx, eeg_notch.chanlocs, ...
    'electrodes','labelpoint','chaninfo', eeg_notch.chaninfo, ...
    'maplimits', [ -5 20 ]) ;
title('Bad channel FASTer',  'FontSize', 20)
colormap(gca,'gray')
%x_lim = get(gca, 'Xlim')
%text(x_lim(1)*1.5,0, bad_chan_label, 'fontweight','bold') %,'FontSize',14)


subplot(3,3,[6,9])
topoplot([],eeg_bad_chan.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', eeg_bad_chan.chaninfo);
title(['Bad channel x INTERP (' num2str(round(length(bad_chan_idx)/eeg_notch.nbchan*100)) '%)'], ...
        'FontSize', 20)


% SAMPLEs - - - - - - -
subplot(6,3,2:3)
imagesc(bad_sample_table)
title(['% Bad Samples = ', num2str(bad_sample_ASRperc)], 'FontSize', 20);
ylabel('iteration');
xlabel('(BAD) SAMPLE  (10k samples = 40 sec)')

bad_sample_ASR2_table = double(~eeg_cleanraw_ASR2.etc.clean_sample_mask);
subplot(6,3,5:6)
imagesc(bad_sample_ASR2_table);
colormap(gca,'gray')
title('BAD segment  after initial ASR repair')

if cfg.do_save_fig
    % fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    %     sprintf('%s_%s_preproc_badSamplesChannels.jpg',subid,dataType));
    fname2save = fullfile(cfg.path.data_prep, ...
        sprintf('%s_preproc_badSamplesChannels',cfg.folder_name));
    print(gcf,fname2save, '-dtiffn','-noui');
end


% ######################################################################
%% INTERPOLATE bad channel identified (with ASR)
% ######################################################################
fprintf('Channel interpolation \n')
eeg_cleanraw_badchan_interp = pop_interp(eeg_cleanraw, ...
    eeg_notch.chanlocs, 'spherical');
    %eeg_struct.chanlocs, 'spherical');
% EEG = eeg_cleanraw_badchan_interp; eeglab redraw


% ######################################################################
%% Average re-referencing
% ######################################################################
fprintf('Average re-referencing \n');
eeg_cleanraw_avgref = pop_reref(eeg_cleanraw_badchan_interp, []);


% ######################################################################
%% Channel pruning before ICA
% ######################################################################
fprintf('Channel pruning before ICA \n')
chan_interp_prune = cfg.prep.chan_interp_prune;
eeg_cleanraw_avgref = pop_select(eeg_cleanraw_avgref, ...
    'nochannel', chan_interp_prune);


% ######################################################################
%% ICA
% ######################################################################
fprintf('Running ICA \n')
eeg_cleanraw_avgref_ICA = pop_runica(eeg_cleanraw_avgref, ...
    'icatype', ...
    'runica', ...
    'extended',1, ...
    'interrupt','on');

%% IClabels for classifying ICA components
fprintf('Identifiy brain and non-brain ICA components with IClabels \n')
eeg_cleanraw_avgref_ICA = pop_iclabel(eeg_cleanraw_avgref_ICA, 'default');
% eeg_cleanraw_avgref_ICA = pop_icflag(eeg_cleanraw_avgref_ICA, ...
%     [0 0.2;0.7 1;0.7 1;NaN NaN;0.7 1;0.7 1;0.7 1]);
eeg_cleanraw_avgref_ICA = pop_icflag(eeg_cleanraw_avgref_ICA, ...
    cfg.prep.iclabels_params);

% find the brain and non-brain classified ICA components
brain_ica_components = find(eeg_cleanraw_avgref_ICA.reject.gcompreject == 0);
non_brain_ica_components = find(eeg_cleanraw_avgref_ICA.reject.gcompreject == 1);
n_icacomponents2reject = length(non_brain_ica_components);
n_icacomponents2keep = length(brain_ica_components);

cfg.brain_ica_components = brain_ica_components;
cfg.non_brain_ica_components = non_brain_ica_components;

%% Plot all ICA components and their probability of being classified as brain
% figure for plotting ICA components and their probability of being classified as brain
close all;
brain_prob_thresh = cfg.prep.brain_prob_thresh;

n_components = size(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications,1);
figure; subplot 121;
set(gcf,'Position',[400 500 1400 800])
%plot(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,1) ,1:n_components);
colormap('jet')
scatter_color = ones([1 n_components]) +10;
scatter_color(1,1:n_components/3) = linspace(1,10,n_components/3);
%scatter(x,y,sz,c,'filled')
scatter(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,1) ,1:n_components, 30, scatter_color,'filled')
hold on; plot([brain_prob_thresh,brain_prob_thresh],ylim,'--');
xlim([0,1]);
grid on;
ylabel('ICA Component');
xlabel('Brain Probability');
title('ICA Component Brain Probability');
view([90 -90])

yLims = ylim;
xLims = xlim;
text(xLims(2)-0.10, yLims(2)-40, sprintf('Brain components = %d',n_icacomponents2keep));
text(xLims(2)-0.15, yLims(2)-40, sprintf('Non-brain components = %d',n_icacomponents2reject));
%text(xLims(2)-0.375, yLims(2)-10,sprintf('Brain components = %d',n_icacomponents2keep));
%text(xLims(2)-0.375, yLims(2)-5,sprintf('Non-brain components = %d',n_icacomponents2reject));

% arrow pointing to EYE-component (3° class)
eye_component_idx = find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,3)*100 > 51);
if ~isempty(eye_component_idx)
    for jj = 1:length(eye_component_idx)
    %     ha = annotation('arrow');  % store the arrow information in ha
    %     ha.Parent = gca;           % associate the arrow the the current axes
    %     % the location in data units
    %     ha.X = [ eye_component_idx(jj), 20 ];    % adjust length and location of arrow 
    %     ha.Y = [ eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(eye_component_idx(jj),3) 0.5 ];      
    % % annotation('textarrow',x_tmp,y_tmp,'String','Eye-IC ??','FontSize',13,'Linewidth',2)
          scatter(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(eye_component_idx(jj),1),eye_component_idx(jj), 80, 'k')
    end
end

% IC classes BAR PLOT:
eye_classes = [ length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,1) > 0.51)), ...
                length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,2) > 0.51)), ...
                length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,3) > 0.51)), ...
                length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,4) > 0.51)), ...
                length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,5) > 0.51)), ...
                length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,6) > 0.51)), ...
                length(find(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,7) > 0.51)) ]
subplot 122
IC_classes = categorical({'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'})
bar(IC_classes, eye_classes)
grid on;
title('IC class prob (> 51%)')

if cfg.do_save_fig
    %fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    %    sprintf('%s_%s_preproc_icacomponent_brain_probabilities.jpg',subid,dataType));
    
    fname2save = fullfile(cfg.path.data_prep, ...
        sprintf('%s_preproc_IC_brain_probabilities',cfg.folder_name));
    print(gcf,fname2save, '-dtiffn','-noui');
end


%% Figure showing all ICA components and their classifications
close all;

% we are interested only to the first 20-30 components
try %if cfg.do_remotePC == 0
    pop_viewprops(eeg_cleanraw_avgref_ICA, 0, [1:28], [2 80], []);
catch
    pop_selectcomps(eeg_cleanraw_avgref_ICA, [1:28] );
end

title('Components 1-28');
if cfg.do_save_fig
    fname2save = fullfile(cfg.path.data_prep, ...
        sprintf('%s_preproc_ica_1_28',cfg.folder_name));
    print(gcf,fname2save, '-dtiffn','-noui');
end


% pop_viewprops(eeg_cleanraw_avgref_ICA, 0, [36:70], [2 80], []);
% title('Components 36-70');
% fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
%     sprintf('%s_%s_preproc_ica_36_70.jpg',subid,dataType));
% print(gcf,fname2save, '-djpeg','-noui');
% 
% pop_viewprops(eeg_cleanraw_avgref_ICA, 0, [71:n_components], [2 80], []);
% title('Components 71-93');
% fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
%     sprintf('%s_%s_preproc_ica_71_%d.jpg',subid,dataType,n_components));
% print(gcf,fname2save, '-djpeg','-noui');


%% Project out all of the non-brain ICA components
fprintf('Projecting out all non-brain ICA components with IClabels \n');
if n_icacomponents2keep~=0
    eeg_cleanraw_avgref_nobadICA = pop_subcomp(eeg_cleanraw_avgref_ICA, non_brain_ica_components);
    %eeg_cleanraw_avgref_nobadICA = pop_subcomp(eeg_cleanraw_avgref_ICA, non_brain_ica_components, 1);
    %eeg_cleanraw_avgref_nobadICA = pop_subcomp_av20210928(eeg_cleanraw_avgref_ICA, non_brain_ica_components, 1);
    
    EEG = eeg_cleanraw_avgref_nobadICA;
else
    fprintf('No ICA components classified as brain! \n');
    fprintf('Final data will be data preprocessed with all steps up until ICA. \n');
    EEG = eeg_cleanraw_avgref_ICA;
end % if n_icacomponents2keep~=0


%% Plot figure showing PSDs before and after ICA denoising
close all;
figure; set(gcf,'color','white');
try
    if n_icacomponents2keep~=0
        subplot(1,3,3);
        pop_spectopo(eeg_cleanraw_avgref_nobadICA, 1, [ ], ...
            'EEG' , 'percent', 50, 'freq', [8 13 20], ...
            'freqrange',[2 80],'electrodes','on');
        sgtitle('PSD Before and After ICA');
    end

        subplot(1,3,2);
        pop_spectopo(eeg_cleanraw_avgref_ICA, 1, [ ], ...
            'EEG' , 'percent', 50, 'freq', [8 13 20], ...
            'freqrange',[2 80],'electrodes','on');

        % PLOT the PSD only for the bad channels
        subplot(1,3,1);
        pop_spectopo(eeg_bad_chan, 1, [ ], ...
            'EEG' , 'percent', 50, 'freq', [8 13 20], ...
            'freqrange',[2 80],'electrodes','on');
        title('PSD of bad (interpolated) channels')

catch error_tmp
    disp(error_tmp)
end
        
if cfg.do_save_fig
    fname2save = fullfile(cfg.path.data_prep, ...
        sprintf('%s_preproc_psd',cfg.folder_name));
    print(gcf,fname2save, '-dtiffn','-noui');
end


%% SAVE
% ===============================================
% KEEP data at 4+1 prep stages: 
% after NOTCH; after ASR repair; after ICA rejection
eeg_prep = [];
eeg_prep.eye_struct = eye_struct;
eeg_prep.eeg_notch = eeg_notch;
eeg_prep.eeg_afterASR_reject = eeg_cleanraw_iter;
eeg_prep.eeg_afterASR_repair = eeg_cleanraw;
eeg_prep.eeg_ICA = eeg_cleanraw_avgref_ICA;
eeg_prep.eeg_nobadICA = eeg_cleanraw_avgref_nobadICA;
eeg_prep.cfg = cfg;
%pop_saveset(eeg_cleanraw_avgref_ICA, 'filename', [ 'eeg_cleanraw_avgref_nobadICA_TODELETE' ]);

cd(cfg.path.data_prep)
fname2save = fullfile(cfg.path.data_prep, ...
         sprintf('%s_eeg_struct_step.mat',cfg.folder_name))
save(fname2save, 'eeg_prep')


% TO INSPECT with VISUAL_ARTIFACTS
% help vis_artifacts
%vis_artifacts(eeg_cleanraw_iter, eeg_notch, 'WindowLength', 400)
%vis_artifacts(eeg_cleanraw, eeg_nobad_chan, 'WindowLength', 400)
%vis_artifacts(eeg_cleanraw_badchan_interp, eeg_notch, 'WindowLength', 400)

%vis_artifacts(eeg_prep.eeg_afterASR, eeg_prep.eeg_notch) %, 'WindowLength', 400)

% !!! NOT WORKING - - - - - - - -
%vis_artifacts(eeg_cleanraw_avgref_ICA, eeg_cleanraw_avgref_nobadICA, 'WindowLength', 120)
% ALTERNATIVE. 
if cfg.do_remotePC == 0
    pop_subcomp_av20210928(eeg_cleanraw_avgref_ICA, non_brain_ica_components, 1);

    if cfg.do_save_fig
        fname2save = fullfile(cfg.path.data_prep, ...
            sprintf('%s_nobadICA_comparison',cfg.folder_name));
        print(gcf,fname2save, '-dtiffn','-noui');
    end
end


%% Save final preproc data to disk
eeg_cleanraw_avgref_nobadICA.cfg = cfg;

cd(fullfile(cfg.path.data_prep))
if n_icacomponents2keep~=0
    %fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    %    sprintf('%s_%s_%s.set',subid, dataType, cfg.pp.final_pp_dn_fstem));
    %pop_saveset(eeg_cleanraw_avgref_nobadICA, 'filename', fname2save);
    pop_saveset(eeg_cleanraw_avgref_nobadICA, 'filename', [ cfg.folder_name '_eeg_cleanraw_avgref_nobadICA' ]);
else
%     fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
%         sprintf('%s_%s_%s.set',subid, dataType, cfg.pp.final_pp_nodn_fstem));
%     pop_saveset(eeg_cleanraw_avgref_ICA, 'filename', fname2save);
    pop_saveset(eeg_cleanraw_avgref_ICA, 'filename', [ cfg.folder_name '_eeg_cleanraw_avgref_ICA' ]);
end

% %% Parameter dump for report
% fprintf('Dumping out all parameters used for preprocessing \n');
% parameter_dump4report(cfg);

fprintf('... Done! \n');

end % function cmi_preproc_land
