%%%
%
%
%
% working with Shay's data


load('arhmm_labels.mat','state_labels');
load('experiment_data_scores.mat','frame_idx','smooth_scores');
load('experiment_data_pca.mat','features');
load('experiment_data_rps.mat','rps');
load('experiment_data_metadata.mat','metadata');
use_session=2;

%%

use_metadata=metadata(use_session);
use_rps=zscore(zscore(rps{use_session})');
use_frame_idx=frame_idx{use_session};
use_labels=state_labels{use_session}(~isnan(use_frame_idx));
use_scores=smooth_scores{use_session}(:,~isnan(use_frame_idx));
load(use_metadata.frames_file,'depth_bounded_rotated');

%%

use_insc_ts=use_metadata.inscopix_frame_idx;

% inscopix formatted?

%%


% changepoints first...


win_size=900;
delta_win=3;
smooth_sig=1;
delta_thresh=.3;
max_lag=90;
bin_smooth=2;
nrands=1e3;

if exist('delta_score.mat','file')~=2
   
        
    delta_score=zeros(size(use_rps,1),length(delta_win));
    %tmp_rps=zscore(zscore(use_rps{i})');
    %tmp_rps=hampel(tmp_rps',5,3)';
    deltas=markolab_deltacoef(use_rps,delta_win); % delta coefficients, lag of 4 frames
    delta_score=sum(abs(deltas)>delta_thresh); % binarize deltas
    h=normpdf([-10:10],0,smooth_sig); % gauss smooth
    delta_score=conv(delta_score,h,'same'); % convolve

    save('delta_score.mat','delta_score');
else
    fprintf('Loading delta scores...\n');
    load('delta_score.mat','delta_score');
end

%%
%

% put the relevant shit on the Inscopix grid

insc_min=min(use_insc_ts(~isnan(use_insc_ts)));
insc_max=max(use_insc_ts(~isnan(use_insc_ts)));

insc_ts=[1:size(neuron_results.C,2)];

insc_scores=nan(size(use_scores,1),numel(insc_ts));
insc_delta=nan(1,numel(insc_ts));
insc_frames=zeros(80,80,numel(insc_ts));

upd=kinect_proctimer(length(insc_ts));

for i =1:length(insc_ts)
    
    % simply average for now
    
    hits=find(use_insc_ts==insc_ts(i));
    if ~isempty(hits)
       insc_scores(:,i)=mean(use_scores(:,hits),2);
       insc_delta(i)=mean(delta_score(hits));
       insc_frames(:,:,i)=mean(depth_bounded_rotated(:,:,hits),3);
    end
    upd(i);
    
end

%%
nanidx=isnan(insc_delta);

for i=1:size(insc_scores,1)
    insc_scores(i,nanidx)=interp1(insc_ts(~nanidx),insc_scores(i,~nanidx),insc_ts(nanidx));
end
insc_delta(nanidx)=interp1(insc_ts(~nanidx),insc_delta(~nanidx),insc_ts(nanidx));

%%
% now run the correlations asshole

% pcs x cells x lags
max_lag=90;
npcs=size(insc_scores,1);
ncells=size(neuron_results.C,1);
lag_mat_pc=zeros(2*max_lag+1,npcs,ncells);
upd=kinect_proctimer(npcs*ncells);
counter=1;
for i=1:npcs
   for j=1:ncells
      lag_mat_pc(:,i,j)=xcorr(zscore(neuron_results.C(j,:)),zscore(insc_scores(i,:)),max_lag,'coeff');
      upd(counter);
      counter=counter+1;
   end
end

%%

lag_mat_delta=zeros(ncells,2*max_lag+1);

for i=1:ncells
  lag_mat_delta(i,:)=xcorr(zscore(neuron_results.C(i,:)),zscore(insc_delta),max_lag,'coeff');
end

%%

% trigger the population on changepoints

[~,changepoints]=findpeaks(insc_delta,'minpeakheight',50);

win_population=zeros(2*max_lag+1,numel(changepoints),ncells);
npoints=numel(insc_delta);
counter=1;

for i=1:length(changepoints)
   left_edge=changepoints(i)-max_lag;
   right_edge=changepoints(i)+max_lag;
   
   if left_edge>0 & right_edge<=npoints
        win_population(:,counter,:)=neuron_results.C(:,left_edge:right_edge)';
        counter=counter+1;
   end
end


% zscore,pca, kmeans (20 looks good!)








