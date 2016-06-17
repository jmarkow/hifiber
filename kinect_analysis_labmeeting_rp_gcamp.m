% load in desired data, from 15781, session_20160513172316

m=matfile('depth_bounded_clean.mat');
use_data=m.depth_bounded_clean; % transfer data to workspace
clear m;

% median filter the data

for i=1:size(use_data,3)
    use_data(:,:,i)=medfilt2(use_data(:,:,i),[9 9]);
end

use_features=double(reshape(use_data,[],size(use_data,3))); % reshape data
movie_fs=30;
ms_per_frame=(1/movie_fs)*1e3;

nidaq_data=kinect_nidaq2mat;
gcamp_trace=nidaq_data(2,:);
gcamp_ts=nidaq_data(end,:);

depth_ts=kinect_read_csv;
depth_ts=depth_ts(:,2);


%%
% johnson-lindenstruass bound

jl_bound_eps=.25;
n_components=round(4*log(size(use_data,3))/(jl_bound_eps^2/2-jl_bound_eps^3/3));
%n_components=400;
rps=kinect_gaussproj(use_features',n_components);
rps=zscore(zscore(rps)');

%%
% get delta coefficients

rps_deltas=markolab_deltacoef(rps,4); % delta coefficients, lag of 4 frames
delta_score=sum(abs(rps_deltas)>.1); % binarize deltas
h=normpdf([-20:20],0,.5); % gauss smooth
delta_score=conv(delta_score,h,'same'); % convolve

thresh=mean(delta_score)+.5*std(delta_score);
idx=1:length(delta_score)-1;

%block_onsets=delta_score(idx)<thresh&delta_score(idx+1)>thresh;
[block_peaks,block_locs]=findpeaks(delta_score,'minpeakheight',thresh);
ms_per_frame=1/movie_fs*1e3;
ibi=diff((block_locs)*ms_per_frame);


bins=[0:100:2e3];
n=histc(ibi,bins);

%% photometry analysis

win_size=2;
new_fs=100;
tau=.2;
nidaq_fs=round(1./mean(diff(gcamp_ts)));

[b,a]=ellip(4,.2,40,[20]/(1e3/2),'low');
gcamp_trace_ds=downsample(filtfilt(b,a,-gcamp_trace),nidaq_fs/new_fs);
gcamp_ts_ds=downsample(gcamp_ts,10);
%gcamp_trace_smooth=markolab_smooth(gcamp_trace_ds(:),round(tau*new_fs),'r','e');
[b2,a2]=butter(3,[1]/(new_fs/2),'low');
gcamp_trace_smooth=filtfilt(b,a,gcamp_trace_ds(:));
%gcamp_trace_smooth=gcamp_trace_ds(:);
gcamp_trace_detrended=fluolab_detrend(gcamp_trace_smooth,'fs',new_fs,'win',win_size);

gcamp_trace_detrended=gcamp_trace_detrended(new_fs*win_size:end);
gcamp_ts_ds=gcamp_ts_ds(new_fs*win_size:end);
gcamp_interp=interp1(gcamp_ts_ds,gcamp_trace_detrended,depth_ts,'linear','extrap');

%% xcorr between the two signals
% bootstrap?

nboots=1e3;
bootvals_max=nan(1,nboots);
bootvals_min=nan(1,nboots);

for i=1:nboots
    scr=markolab_phase_scramble_1d(zscore(gcamp_interp),0);
    r=xcorr(zscore(scr),zscore(delta_score),25,'coeff');
    bootvals_max(i)=max(r);
    bootvals_min(i)=min(r);
end

[obs_r,lags]=xcorr(zscore(gcamp_interp),zscore(delta_score),25,'coeff');

%%
%gcamp_score=zscore(gcamp_trace_detrended);
%gcamp_score(gcamp_score<0)=0;
%gcamp_score=gcamp_score;
[pk.vals,pk.locs]=findpeaks(gcamp_trace_detrended,'minpeakheight',1,'minpeakdist',5);
gcamp_binary=zeros(size(gcamp_trace_detrended));
gcamp_binary(pk.locs)=1;

% compare w/ block locations

% for each calcium peak, bin delta score (or the other way around)

% get peaks, then rise times?


%%
figure();
markolab_stairplot(n,bins,'method','a');

%%
figure();
ax(1)=subplot(3,1,1);
imagesc(depth_ts-min(depth_ts),[],filter(ones(3,1)/3,1,rps')');
colormap(bone)
axis off;
ax(2)=subplot(3,1,2);
plot(depth_ts-min(depth_ts),conv(delta_score,normpdf([-10:10],0,1.5),'same'));
set(ax(2),'ydir','rev');
axis off;
ax(3)=subplot(3,1,3);
plot(gcamp_ts_ds-min(depth_ts),gcamp_trace_detrended);
ylim([-1 6]);
box off;
set(gca,'TickDir','out','TickLength',[.025 .025]);
linkaxes(ax,'x');
% decent example from 15781 20160513172316
xlim([11 34]);


%%
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
