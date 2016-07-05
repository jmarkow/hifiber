%%
% analysis of 15781, D1-cre

m=matfile('depth_bounded_rotated.mat');
use_data=m.depth_bounded_rotated; % transfer data to workspace
use_mask=m.depth_bounded_cable_mask_rotated;
rankcut=20;
nboots=1e3;

% get timestamps and remaining data

corr_limits=[-1 1 0 .3];
movie_fs=30;
ms_per_frame=(1/movie_fs);

nidaq_data=kinect_nidaq2mat;
photometry_traces=nidaq_data(1:end-1,:);
photometry_ts=nidaq_data(end,:);

depth_ts=kinect_read_csv;
depth_ts=depth_ts(:,2);

%%

features=single(reshape(use_data,[],size(use_data,3)));
mask=reshape(use_mask,[],size(use_mask,3));
mask2=features<15;
features(~mask)=nan;
features(mask2)=0;

% clean up behavioral data


[U S V mu]=kinect_randsvd_power_iter(features,[],[],3);
recon=U(:,1:rankcut)*S(1:rankcut,1:rankcut)*V(:,1:rankcut)'+mu;
[rps,delta_score,changepoints,changethresh]=kinect_analysis_changepoints(recon','delta_thresh',.05,'smooth_sig',2);
recon=reshape(recon,sqrt(size(recon,1)),sqrt(size(recon,1)),[]);

%%
% process photometry

[gcamp_data,gcamp_ts]=kinect_analysis_proc_photometry(photometry_traces(1,:)+photometry_traces(2,:)*1e2,photometry_ts,'dff',1,'smooth_tau',.2);
gcamp_interp=interp1(gcamp_ts,gcamp_data,depth_ts,'linear','extrap');
[rcamp_data,rcamp_ts]=kinect_analysis_proc_photometry(photometry_traces(3,:)+photometry_traces(4,:)*1e2,photometry_ts,'dff',1,'smooth_tau',.2);
rcamp_interp=interp1(rcamp_ts,rcamp_data,depth_ts,'linear','extrap');

% xcorr between the two signals

[bootvals_gcamp.max,bootvals_gcamp.min,obs_r_gcamp,lags]=kinect_analysis_xcorr_bootstr(zscore(gcamp_interp),zscore(delta_score));
[bootvals_rcamp.max,bootvals_rcamp.min,obs_r_rcamp]=kinect_analysis_xcorr_bootstr(zscore(rcamp_interp),zscore(delta_score));


%%

figs.example=figure();
ax(1)=subplot(4,1,1);
imagesc(depth_ts-min(depth_ts),[],squeeze(mean(recon(40:60,:,:))));
caxis([0 30]);
axis off;
ax(2)=subplot(4,1,2);
imagesc(depth_ts-min(depth_ts),[],filter(ones(5,1)/5,1,rps')');
colormap(bone)
axis off;
ax(3)=subplot(4,1,3);
plot(depth_ts-min(depth_ts),conv(delta_score,normpdf([-10:10],0,1.5),'same'));
set(ax(3),'ydir','rev');
ylim([changethresh max(delta_score)+eps]);
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
figs.xcorr=figure();
subplot(1,2,1);
ydata=[repmat(prctile(bootvals_gcamp.max,99),size(lags));repmat(prctile(bootvals_gcamp.min,1),size(lags))];
h1=markolab_shadeplot(lags*ms_per_frame,ydata,[.7 .7 .7],'none',2);
hold on;
plot(lags*ms_per_frame,obs_r_gcamp);
plot(lags*ms_per_frame,repmat(prctile(bootvals_gcamp.max,99),size(lags)),'k--');
plot([0 0],[-1 1],'k-');

set(gca,'Layer','Top','xtick',[corr_limits(1):.5:corr_limits(2)],'YTick',[0:.1:corr_limits(4)]);
axis([corr_limits]);
uistack(h1,'bottom');


subplot(1,2,2);
ydata=[repmat(prctile(bootvals_rcamp.max,99),size(lags));repmat(prctile(bootvals_rcamp.min,1),size(lags))];
h1=markolab_shadeplot(lags*ms_per_frame,ydata,[.7 .7 .7],'none',2);
hold on;
plot(lags*ms_per_frame,obs_r_rcamp);
plot(lags*ms_per_frame,repmat(prctile(bootvals_rcamp.max,99),size(lags)),'k--');
plot([0 0],[-1 1],'k-');

set(gca,'Layer','Top','xtick',[corr_limits(1):.5:corr_limits(2)],'YTick',[0:.1:corr_limits(4)]);
axis([corr_limits]);
uistack(h1,'bottom');
set(figs.xcorr,'position',[300 300 350 150]);