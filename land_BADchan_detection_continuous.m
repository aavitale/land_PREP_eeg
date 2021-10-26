function [ FASTbadChans, chan_properties, z_score ] = land_BADchan_detection_continuous(EEG, z_threshold);

%% STEP 7 of FASTER -> to find bad channels

    % First check whether reference channel (i.e. zeroed channels) is present in data
    % reference channel is needed to run faster
    ref_chan=[]; FASTbadChans=[]; all_chan_bad_FAST=0;
    ref_chan=find(any(EEG.data, 2)==0);
    if numel(ref_chan)>1
        error(['There are more than 1 zeroed channel (i.e. zero value throughout recording) in data.'...
            ' Only reference channel should be zeroed channel. Delete the zeroed channel/s which is not reference channel.']);
    elseif numel(ref_chan)==1
        
        % array with: 
        % 1°column: Mean CORRELATION between each channel and all other channels
        % 2°column: VARIANCE of the channel
        % 3°coumn: HURST exponent        
        
        eeg_chans = 1:EEG.nbchan;
        chan_properties = channel_properties(EEG, eeg_chans, ref_chan); % run faster

        % compute the Z-score 
        % and identify channels with +/- 3std in one of those metrics
        %[ FASTbadIdx, z_score ] = min_z(chan_properties);
        
        if ~exist('z_threshold')
            z_threshold = 3 %2.5 %3
        end
        rejection_options.measure=logical(ones(1,size(chan_properties,2)));
        z_score=chan_properties-repmat(mean(chan_properties,1),size(chan_properties,1),1);
        z_score=z_score./repmat(std(z_score,[],1),size(chan_properties,1),1);
        z_score(isnan(z_score))=0;
        all_l = abs(z_score) > repmat(z_threshold, size(chan_properties,1),1);
        FASTbadIdx = sum(all_l,2) %sumany(all_l(:,rejection_options.measure),2);
        
        %FASTbadChans=find(FASTbadIdx==1);
        FASTbadChans=find(FASTbadIdx>0);
        FASTbadChans=FASTbadChans(FASTbadChans~=ref_chan);
        
        %reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
        
        EEG = eeg_checkset(EEG);
        channels_analysed=EEG.chanlocs; % keep full channel locations to use later for interpolation of bad channels
%     elseif numel(ref_chan)==0
%         warning('Reference channel is not present in data. Cz channel will be used as reference channel.');
%         ref_chan=find(strcmp({EEG.chanlocs.labels}, 'Cz')); % find Cz channel index
%         EEG_copy=[];
%         EEG_copy=EEG; % make a copy of the dataset
%         EEG_copy = pop_reref( EEG_copy, ref_chan,'keepref','on'); % rerefer to Cz in copied dataset
%         EEG_copy = eeg_checkset(EEG_copy);
%         list_properties = channel_properties(EEG_copy, 1:EEG_copy.nbchan, ref_chan); % run faster on copied dataset
%         FASTbadIdx=min_z(list_properties);
%         FASTbadChans=find(FASTbadIdx==1);
%         channels_analysed=EEG.chanlocs;
%         reference_used_for_faster{subject}={EEG.chanlocs(ref_chan).labels};
    end
    
    
% figure; hold on
% badchan_idx = [ 1 2 3 ]
% for jj = 1:length(badchan_idx)
%     subplot(2,2,jj);
%     pop_prop( EEG, 1, [badchan_idx(jj)], NaN, {'freqrange',[2 45] });
% end


%% sub-FUNCTIONs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function list_properties = channel_properties(EEG,eeg_chans,ref_chan)

    if ~isempty(ref_chan) && length(ref_chan)==1
        pol_dist=distancematrix(EEG,eeg_chans);
        [s_pol_dist dist_inds] = sort(pol_dist(ref_chan,eeg_chans));
        [s_inds idist_inds] = sort(dist_inds);
    end

    % TEMPORAL PROPERTIES
    % 1 Mean correlation between each channel and all other channels

    % Ignore zeroed channels (ie reference channels) to avoid NaN problems
    ignore = [];
    datacorr = EEG.data;
    for u = eeg_chans
        if max(EEG.data(u,:))==0 && min(EEG.data(u,:))==0
            ignore=[ignore u];
        end
    end
    
    measure = 1;

    % Calculate correlations
    calc_indices=setdiff(eeg_chans,ignore);
    ignore_indices=intersect(eeg_chans,ignore);
    corrs = abs(corrcoef(EEG.data(setdiff(eeg_chans,ignore),:)'));
    mcorrs=zeros(size(eeg_chans));
    for u=1:length(calc_indices)
        mcorrs(calc_indices(u))=mean(corrs(u,:));
    end
    mcorrs(ignore_indices)=mean(mcorrs(calc_indices));

    % Quadratic correction for distance from reference electrode

    if (~isempty(ref_chan) && length(ref_chan)==1)
        p = polyfit(s_pol_dist,mcorrs(dist_inds),2);
        fitcurve = polyval(p,s_pol_dist);
        corrected = mcorrs(dist_inds) - fitcurve(idist_inds);

        list_properties(:,measure) = corrected;
    else
        list_properties(:,measure) = mcorrs(dist_inds);
    end
    
    
    % 2 VARIANCE of the channels
    measure = measure + 1;
    
    vars = var(EEG.data(eeg_chans,:)');
    vars(~isfinite(vars))=mean(vars(isfinite(vars)));
    % Quadratic correction for distance from reference electrode

    if (~isempty(ref_chan) && length(ref_chan)==1)
        p = polyfit(s_pol_dist,vars(dist_inds),2);
        fitcurve = polyval(p,s_pol_dist);
        corrected = vars - fitcurve(idist_inds);

        list_properties(:,measure) = corrected;
    else
        list_properties(:,measure) = vars;
    end
    measure = measure + 1;

    
    % 3 Hurst exponent
    for u=1:length(eeg_chans)
        list_properties(u,measure) = hurst_exponent(EEG.data(eeg_chans(u),:));
    end

    for u = 1:size(list_properties,2)
        list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
        list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
    end
end


function [distmatrixpol distmatrixxyz distmatrixproj] = distancematrix(EEG,eeg_chans)

    num_chans = size(EEG.data,1);
    distmatrix = zeros(length(eeg_chans),length(eeg_chans));
    distmatrixpol = [];
    for chan2tst = eeg_chans;
        for q=eeg_chans
            distmatrixpol(chan2tst,q)=sqrt(((EEG.chanlocs(chan2tst).radius^2)+(EEG.chanlocs(q).radius^2))-(2*((EEG.chanlocs(chan2tst).radius)*...
                (EEG.chanlocs(q).radius)*cosd(EEG.chanlocs(chan2tst).theta - EEG.chanlocs(q).theta))));%calculates the distance between electrodes using polar format
        end
    end

    locs = EEG.chanlocs;
    for u = eeg_chans
        if ~isempty(locs(u).X)
            Xs(u) = locs(u).X;
        else
            Xs(u) = 0;
        end
        if ~isempty(locs(u).Y)
            Ys(u) = locs(u).Y;
        else
            Ys(u) = 0;
            end
        if ~isempty(locs(u).Z)
            Zs(u) = locs(u).Z;
        else
            Zs(u) = 0;
        end
    end
    Xs = round2(Xs,6);
    Ys = round2(Ys,6);
    Zs = round2(Zs,6);

    for u = eeg_chans
        for v=eeg_chans
            distmatrixxyz(u,v) = dist(Xs(u),Xs(v))+dist(Ys(u),Ys(v))+dist(Zs(u),Zs(v));
        end
    end
    D = max(max(distmatrixxyz));
    distmatrixproj = (pi-2*(acos(distmatrixxyz./D))).*(D./2);
        function d = dist(in1,in2)
            d = sqrt(abs(in1.^2 - in2.^2));
        end

        function num = round2(num,decimal)
            num = num .* 10^decimal;
            num = round(num);
            num = num ./ 10^decimal;
        end
end


% The Hurst exponent
%--------------------------------------------------------------------------
% This function does dispersional analysis on a data series, then does a 
% Matlab polyfit to a log-log plot to estimate the Hurst exponent of the 
% series.
%
% This algorithm is far faster than a full-blown implementation of Hurst's
% algorithm.  I got the idea from a 2000 PhD dissertation by Hendrik J 
% Blok, and I make no guarantees whatsoever about the rigor of this approach
% or the accuracy of results.  Use it at your own risk.
%
% Bill Davidson
% 21 Oct 2003

function [hurst] = hurst_exponent(data0)   % data set

    data=data0;         % make a local copy

    [M,npoints]=size(data0);

    yvals=zeros(1,npoints);
    xvals=zeros(1,npoints);
    data2=zeros(1,npoints);

    index=0;
    binsize=1;

    while npoints>4

        y=std(data);
        index=index+1;
        xvals(index)=binsize;
        yvals(index)=binsize*y;

        npoints=fix(npoints/2);
        binsize=binsize*2;
        for ipoints=1:npoints % average adjacent points in pairs
            data2(ipoints)=(data(2*ipoints)+data((2*ipoints)-1))*0.5;
        end
        data=data2(1:npoints);

    end % while

    xvals=xvals(1:index);
    yvals=yvals(1:index);

    logx=log(xvals);
    logy=log(yvals);

    p2=polyfit(logx,logy,1);
    hurst=p2(1); % Hurst exponent is the slope of the linear fit of log-log plot
end


end


