%%
% analysis of 15781, D1-cre

m=matfile('depth_bounded_rotated.mat');
use_data=m.depth_bounded_rotated; % transfer data to workspace
use_mask=m.depth_bounded_cable_mask_rotated;
rankcut=20;
nboots=1e3;

% get timestamps and remaining data

movie_fs=30;
ms_per_frame=(1/movie_fs)*1e3;

nidaq_data=kinect_nidaq2mat;
photometry_traces=nidaq_data(1:end-1,:);
photometry_ts=nidaq_data(end,:);

depth_ts=kinect_read_csv;
depth_ts=depth_ts(:,2);

features=single(reshape(use_data,[],size(use_data,3)));
mask=reshape(use_mask,[],size(use_mask,3));
mask2=features<15;
features(~mask)=nan;
features(mask2)=0;

%% clean up behavioral data


[U S V mu]=kinect_randsvd_power_iter(features,[],[],3);

%%


recon=U(:,1:rankcut)*S(1:rankcut,1:rankcut)*V(:,1:rankcut)'+mu;
[rps,delta_score,changepoints,changethresh]=kinect_analysis_changepoints(recon','delta_thresh',.2,'smooth_sig',2);
recon=reshape(recon,sqrt(size(recon,1)),sqrt(size(recon,1)),[]);

%%
% process photometry

[gcamp_data,gcamp_ts]=kinect_analysis_proc_photometry(photometry_traces(1,:),photometry_ts,'dff',0,'smooth_tau',.4,'detrend_win_size',3);
[rcamp_data,gcamp_ts]=kinect_analysis_proc_photometry(photometry_traces(3,:),photometry_ts,'dff',0,'smooth_tau',.4,'detrend_win_size',3);

%gcamp_data=gcamp_data/1e2;
% to_del=find(depth_ts<min(gcamp_ts));
% depth_ts(to_del)=[];
% rps(:,to_del)=[];
% delta_score(to_del)=[];
% recon(:,:,to_del)=[];

gcamp_interp=interp1(gcamp_ts,gcamp_data,depth_ts,'linear','extrap');
rcamp_interp=interp1(gcamp_ts,rcamp_data,depth_ts,'linear','extrap');

%%

% get gcamp peaks

gcamp_score=zscore(gcamp_data);
% gcamp_score(gcamp_score<0)=0;
% gcamp_score=gcamp_score.^2;
idx=[1:length(gcamp_score)-1];
gcamp_idx=[0;gcamp_score(idx)<.5&gcamp_score(idx+1)>.5];
gcamp_locs=find(gcamp_idx);

rcamp_wins_gcamp_pks=kinect_analysis_win_data(gcamp_locs,150,zscore(rcamp_data));

% bootstrap the analysis with random sampling 

rcamp_bootstr_wins=zeros(size(rcamp_wins_gcamp_pks,1),nboots);
pool=101:length(gcamp_data)-101;

for i=1:nboots
    tmp=kinect_analysis_win_data(randsample(pool,length(gcamp_locs)),150,zscore(rcamp_data));
    rcamp_bootstr_wins(:,i)=mean(tmp,2);
end

rcamp_score=zscore(rcamp_data);
% gcamp_score(gcamp_score<0)=0;
% gcamp_score=gcamp_score.^2;
rcamp_idx=[0;rcamp_score(idx)<.5&rcamp_score(idx+1)>.5];
rcamp_locs=find(rcamp_idx);

gcamp_wins_rcamp_pks=kinect_analysis_win_data(rcamp_locs,150,zscore(gcamp_data));

% bootstrap the analysis with random sampling 

gcamp_bootstr_wins=zeros(size(gcamp_wins_rcamp_pks,1),nboots);

for i=1:nboots
    tmp=kinect_analysis_win_data(randsample(pool,length(rcamp_locs)),150,zscore(gcamp_data));
    gcamp_bootstr_wins(:,i)=mean(tmp,2);
end

%% xcorr between the two signals
% bootstrap?

bootvals_max_gcamp=nan(1,nboots);
bootvals_min_gcamp=nan(1,nboots);

for i=1:nboots
    scr=markolab_phase_scramble_1d(zscore(gcamp_interp),0);
    r=xcorr(zscore(scr),zscore(delta_score),100,'coeff');
    bootvals_max_gcamp(i)=max(r);
    bootvals_min_gcamp(i)=min(r);
end

[gcamp_obs_r,lags]=xcorr(zscore(gcamp_interp),zscore(delta_score),100,'coeff');

bootvals_max_rcamp=nan(1,nboots);
bootvals_min_rcamp=nan(1,nboots);

for i=1:nboots
    scr=markolab_phase_scramble_1d(zscore(rcamp_interp),0);
    r=xcorr(zscore(scr),zscore(delta_score),100,'coeff');
    bootvals_max_rcamp(i)=max(r);
    bootvals_min_rcamp(i)=min(r);
end

[rcamp_obs_r,lags]=xcorr(zscore(rcamp_interp),zscore(delta_score),100,'coeff');


% and check correlations with each other


%%

% get peaks in delta score to trigger calcium activity
% 
% idx=[1:length(delta_score)-1];
% beh_idx=[0 delta_score(idx)<changethresh*.5&delta_score(idx+1)>changethresh*.5];
% beh_locs=find(beh_idx);
% 
% % get photometry index for each beh_loc
% 
% phot_locs=nan(size(beh_locs));
% 
% for i=1:length(beh_locs)
%     
%     [val,tmp]=min(abs(depth_ts(beh_locs(i))-gcamp_ts));
%     if val<.05
%         phot_locs(i)=tmp;
%     end
%     
% end


%%

figure();
ax(1)=subplot(4,1,1);
imagesc(depth_ts-min(depth_ts),[],squeeze(mean(recon(40:60,:,:))));
caxis([0 30]);
axis off;
ax(2)=subplot(4,1,2);
imagesc(depth_ts-min(depth_ts),[],filter(ones(3,1)/3,1,rps')');
colormap(bone)
axis off;
ax(3)=subplot(4,1,3);
plot(depth_ts-min(depth_ts),conv(delta_score,normpdf([-10:10],0,1.5),'same'));
set(ax(3),'ydir','rev');
ylim([changethresh max(delta_score)+eps]);
axis off;
ax(4)=subplot(4,1,4);
gcamp_data=gcamp_data-min(gcamp_data);
rcamp_data=rcamp_data-min(rcamp_data);
plot(gcamp_ts-min(depth_ts),gcamp_data./max(gcamp_data));
hold on;
plot(gcamp_ts-min(depth_ts),rcamp_data./max(rcamp_data));
ylim([-.01 1]);
box off;
set(gca,'TickDir','out','TickLength',[.025 .025]);
linkaxes(ax,'x');
% decent example from 16639 20160618193331
xlim([225 250]);


%%

%
figure();
ydata=[repmat(prctile(bootvals_max_gcamp,99),size(lags));repmat(prctile(bootvals_min_gcamp,1),size(lags))];
h1=markolab_shadeplot(lags*ms_per_frame,ydata,[.7 .7 .7],'none',2);
hold on;
plot(lags*ms_per_frame,gcamp_obs_r);
plot(lags*ms_per_frame,repmat(prctile(bootvals_max_gcamp,99),size(lags)),'k--');
plot([0 0],[-.2 .2],'k-');

set(gca,'Layer','Top','xtick',[-400:100:400],'YTick',[0:.02:.2]);
axis([-400 400 0 .2]);

xticks=repmat(get(gca,'XTick'),[2 1]);
yticks=repmat(get(gca,'YTick'),[2 1]);
ylimits=ylim();
xlimits=xlim();
h2=plot(xticks,repmat([ylimits(1);ylimits(2)],[1 size(xticks,2)]),'k-','color',[.8 .8 .8]);
h3=plot(repmat([xlimits(1);xlimits(2)],[1 size(yticks,2)]),yticks,'k-','color',[.8 .8 .8]);

uistack(h2,'bottom');
uistack(h3,'bottom');
uistack(h1,'bottom');


%%

figure();
rcamp_gcamp_pks_lags=[-150:150];
ydata=[repmat(prctile(max(rcamp_bootstr_wins),99),size(rcamp_gcamp_pks_lags));repmat(prctile(min(rcamp_bootstr_wins),1),size(rcamp_gcamp_pks_lags))];
h1=markolab_shadeplot(rcamp_gcamp_pks_lags*ms_per_frame,ydata,[.7 .7 .7],'none',2);
hold on;
plot(rcamp_gcamp_pks_lags*ms_per_frame,mean(rcamp_wins_gcamp_pks,2));
plot([0 0],[-.6 .3],'k-');

figure();
rcamp_gcamp_pks_lags=[-150:150];
ydata=[repmat(prctile(max(gcamp_bootstr_wins),99),size(rcamp_gcamp_pks_lags));repmat(prctile(min(gcamp_bootstr_wins),1),size(rcamp_gcamp_pks_lags))];
h1=markolab_shadeplot(rcamp_gcamp_pks_lags*ms_per_frame,ydata,[.7 .7 .7],'none',2);
hold on;
plot(rcamp_gcamp_pks_lags*ms_per_frame,mean(gcamp_wins_rcamp_pks,2));
plot([0 0],[-.6 .3],'k-');
