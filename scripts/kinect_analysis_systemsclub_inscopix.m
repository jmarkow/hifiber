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

camera_fs=30;

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

[~,changepoints]=findpeaks(insc_delta,'minpeakheight',125);

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
%%
changemat=squeeze(mean(win_population(80:110,:,:)));
[coeff score]=pca(zscore(changemat));

% do this the old-fashioned way

clust_choice=[5:30];
bic=[];
upd=kinect_proctimer(length(clust_choice));
options=statset('MaxIter',1e3);

for i=1:length(clust_choice);
    guess=kmeans(zscore(score(:,1:10)),clust_choice(i),'replicates',5);
    model_obj{i}=fitgmdist(zscore(score(:,1:10)),clust_choice(i),'start',guess,'regularization',1e-5,'options',options);
    bic(i)=model_obj{i}.BIC;
    upd(i);
end

save('activity_cluster.mat','model_obj','bic','clust_choice');

%%

[~,loc]=min(bic);
use_model=model_obj{loc};
cidx=cluster(use_model,zscore(score(:,1:10)));
[~,cplotidx]=sort(cidx,'descend');

%%

figs.all_changepoint_activity=figure('PaperPositionMode','Auto','Position',[100 100 700 500]);
whitebg(figs.all_changepoint_activity);
imagesc(zscore(changemat)');
axis off
caxis([0 10]);
colormap(hot);
set(figs.all_changepoint_activity,'color',[0 0 0],'InvertHardcopy','off');
ylimits=ylim()
markolab_multi_fig_save(figs.all_changepoint_activity,'~/Desktop/quickfigs','systemsclub_inscopix_allchangepointactivity','eps,png,fig','renderer','painters');

figs.all_changepoint_activity_cov=figure('PaperPositionMode','Auto','Position',[100 100 700 500]);
whitebg(figs.all_changepoint_activity);
imagesc(cov(zscore(changemat)'));
axis off
caxis([0 .5]);
colormap(hot);
set(figs.all_changepoint_activity_cov,'color',[0 0 0],'InvertHardcopy','off');
ylimits=ylim()
markolab_multi_fig_save(figs.all_changepoint_activity_cov,'~/Desktop/quickfigs','systemsclub_inscopix_allchangepointactivity_cov','eps,png,fig','renderer','painters');

%%

figs.all_changepoint_activity_sorted=figure('PaperPositionMode','Auto','Position',[100 100 700 500]);
whitebg(figs.all_changepoint_activity);
imagesc(imgaussfilt(zscore(changemat(cplotidx,:))',[.5 5]));
axis off
caxis([0 5]);
colormap(hot);
set(figs.all_changepoint_activity,'color',[0 0 0],'InvertHardcopy','off');
ylimits=ylim()
markolab_multi_fig_save(figs.all_changepoint_activity_sorted,'~/Desktop/quickfigs','systemsclub_inscopix_allchangepointactivity_sorted','eps,png,fig','renderer','painters');

figs.all_changepoint_activity_sorted_cov=figure('PaperPositionMode','Auto','Position',[100 100 700 500]);
whitebg(figs.all_changepoint_activity);
imagesc(cov(zscore(changemat(cplotidx,:))'));
axis off
caxis([0 .5]);
colormap(hot);
set(figs.all_changepoint_activity_cov,'color',[0 0 0],'InvertHardcopy','off');
ylimits=ylim()
markolab_multi_fig_save(figs.all_changepoint_activity_sorted_cov,'~/Desktop/quickfigs','systemsclub_inscopix_allchangepointactivity_sorted_cov','eps,png,fig','renderer','painters');


%%
% first show all activity

figs.all_activity=figure('PaperPositionMode','Auto','Position',[100 100 700 500]);
imagesc(zscore(neuron_results.C')');
caxis([0 10]);
hold on;
ylimits=ylim();
plot([changepoints;changepoints],[ones(size(changepoints))*ylimits(1);ones(size(changepoints))*ylimits(2)],'w-');
set(figs.all_activity,'InvertHardcopy','off');
markolab_multi_fig_save(figs.all_activity,'~/Desktop/quickfigs','systemsclub_inscopix_allpluschangepoints','eps,png,fig','renderer','painters');
    

%%
figs.changepoints_corr=figure('PaperPositionMode','Auto','Position',[100 100 300 250]);
whitebg(figs.changepoints_corr)
thresh=std(neuron_results.C,[],2);
bin_activity=neuron_results.C>repmat(thresh,[1 size(neuron_results.C,2)]);
[r,lags]=xcorr(zscore(sum(bin_activity)),zscore(insc_delta),max_lag,'coeff');

%[r,lags]=
g=gramm('x',lags/camera_fs,'y',r);
g.geom_line();
g.axe_property('ylim',[-.08 .08],'ytick',[-.08 0 .08],'xtick',[-3 0 3],'xlim',[-3 3],'TickDIr','out','TickLength',[.015 .015]);

g.draw();
set(figs.changepoints_corr,'color',[0 0 0],'InvertHardCopy','off');
offsetAxes(g.facet_axes_handles);
markolab_multi_fig_save(figs.changepoints_corr,'~/Desktop/quickfigs','systemsclub_inscopix_changepointscorr','eps,png,fig','renderer','painters');


%%
% write out two movies, make sure they are sync'ed up

tf_listing=dir('recording*.tif');

% last is first unfortunately
nframes=1e3;
% v=VideoWriter('rawdata.mp4','MPEG-4')
% v.FrameRate=30;
% v.Quality=100;

frames=zeros(1080,1440,1e3,'int16');
upd=kinect_proctimer(nframes);
for i=1:nframes
    frames(:,:,i)=imread(tf_listing(end).name,i);
    upd(i);
end




