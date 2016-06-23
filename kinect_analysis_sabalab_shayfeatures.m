m=matfile('depth_bounded_rotated.mat');
m2=matfile('depth_masked.mat');
use_data=m.depth_bounded_rotated; % transfer data to workspace
use_mask=m.depth_bounded_cable_mask_rotated;
frame_idx=m2.frame_idx;
rankcut=20;
nboots=1e3;
movie_fs=30;
ms_per_frame=(1/movie_fs)*1e3;


nidaq_data=kinect_nidaq2mat;
depth_ts=kinect_read_csv;
depth_ts=depth_ts(:,2);

depth_ts=depth_ts(frame_idx);

%%

% get flips

fid=fopen('flips.txt');
flips=fscanf(fid,'%f',[1 inf]);
fclose(fid);

for i=1:length(flips)
   use_data(:,:,flips(i):end)=fliplr(use_data(:,:,flips(i):end)); 
   use_mask(:,:,flips(i):end)=fliplr(use_mask(:,:,flips(i):end));
end

%%

% Inscopix frame onsets (NiDAQ clock)

idx=1:size(nidaq_data,2)-1;
insc_frame_onsets=nidaq_data(1,idx)<2.5&nidaq_data(1,idx+1)>2.5;
insc_frame_onsets_nidaq=nidaq_data(2,insc_frame_onsets);

% now get Kinect timestamps from nidaq and match up

insc_frame_idx=zeros(1,length(depth_ts));
dist=zeros(1,length(depth_ts));
for i=1:length(depth_ts)
    [dist(i),insc_frame_idx(i)]=min(abs(depth_ts(i)-insc_frame_onsets_nidaq));
end
insc_frame_idx(abs(dist)>.018)=nan;
% get timestamps and remaining data

%%

features=single(reshape(use_data,[],size(use_data,3)));
mask=reshape(use_mask,[],size(use_mask,3));
mask2=features<15;
features(~mask)=nan;
features(mask2)=0;

%% clean up behavioral data

[U S V mu]=kinect_randsvd_power_iter(features,[],[],5);

%%

recon=U(:,1:rankcut)*S(1:rankcut,1:rankcut)*V(:,1:rankcut)'+mu;
[rps,delta_score,changepoints,changethresh]=kinect_analysis_changepoints(recon','delta_thresh',.15,'smooth_sig',2);
recon=reshape(recon,sqrt(size(recon,1)),sqrt(size(recon,1)),[]);

%%
[Uraw Sraw Vraw]=kinect_randsvd_power(single(reshape(use_data,[],size(use_data,3))),400,2);

figure();
subplot(1,2,1);
kinect_eigenmontage(U,'k',20,'clip',75);
subplot(1,2,2);
kinect_eigenmontage(Uraw,'k',20,'clip',75);
% get usv from raw data

%%

figure();
ax(1)=subplot(3,1,1);
imagesc(depth_ts-min(depth_ts),[],squeeze(mean(recon(40:60,:,:))));
caxis([0 30]);
axis off;
ax(2)=subplot(3,1,2);
imagesc(depth_ts-min(depth_ts),[],filter(ones(5,1)/5,1,rps')');
colormap(bone)
axis off;
ax(3)=subplot(3,1,3);
plot(depth_ts-min(depth_ts),conv(delta_score,normpdf([-10:10],0,1.5),'same'));
set(ax(3),'ydir','rev');
ylim([changethresh max(delta_score)+eps]);
axis off;
linkaxes(ax,'x');
% decent example from 15781 20160513172316
xlim([70 92]);