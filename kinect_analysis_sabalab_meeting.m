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


[U S V mu]=kinect_randsvd_power_iter(features);
recon=U(:,1:rankcut)*S(1:rankcut,1:rankcut)*V(:,1:rankcut)'+mu;
[rps,delta_score,changepoints]=kinect_analysis_changepoints(recon');
recon=reshape(recon,sqrt(size(recon,1)),sqrt(size(recon,1)),[]);
% process photometry

[gcamp_data,gcamp_ts]=kinect_analysis_proc_photometry(photometry_traces(1,:),photometry_ts,'dff',0);
%gcamp_data=gcamp_data/1e2;
gcamp_interp=interp1(gcamp_ts,gcamp_data,depth_ts,'linear','extrap');

% then plot everything


%% xcorr between the two signals
% bootstrap?

bootvals_max=nan(1,nboots);
bootvals_min=nan(1,nboots);

for i=1:nboots
    scr=markolab_phase_scramble_1d(zscore(gcamp_interp),0);
    r=xcorr(zscore(scr),zscore(delta_score),100,'coeff');
    bootvals_max(i)=max(r);
    bootvals_min(i)=min(r);
end

[obs_r,lags]=xcorr(zscore(gcamp_interp),zscore(delta_score),100,'coeff');


%%

figure();
ax(1)=subplot(4,1,1);
imagesc(depth_ts-min(depth_ts),[],squeeze(mean(recon(40:60,:,:))));
colormap(winter);
axis off;
ax(2)=subplot(4,1,2);
imagesc(depth_ts-min(depth_ts),[],filter(ones(3,1)/3,1,rps')');
colormap(bone)
axis off;
ax(3)=subplot(4,1,3);
plot(depth_ts-min(depth_ts),conv(delta_score,normpdf([-10:10],0,1.5),'same'));
set(ax(3),'ydir','rev');
axis off;
ax(4)=subplot(4,1,4);
plot(gcamp_ts-min(depth_ts),gcamp_data);
ylim([-1 6]);
box off;
set(gca,'TickDir','out','TickLength',[.025 .025]);
linkaxes(ax,'x');
% decent example from 15781 20160513172316
xlim([11 34]);


%%

%
figure();
ydata=[repmat(prctile(bootvals_max,99),size(lags));repmat(prctile(bootvals_min,1),size(lags))];
h1=markolab_shadeplot(lags*ms_per_frame,ydata,[.7 .7 .7],'none',2);
hold on;
plot(lags*ms_per_frame,obs_r);
plot(lags*ms_per_frame,repmat(prctile(bootvals_max,99),size(lags)),'k--');
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
