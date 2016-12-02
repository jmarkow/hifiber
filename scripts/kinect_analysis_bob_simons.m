%%%

% take some shit, do some shit


% assume pcs are loaded in, use for reconstruction


% use day

use_session=2;

% reconstruct the frames?

load(metadata(use_session).frames_file,'depth_bounded_rotated','depth_bounded_cable_mask_rotated');
orientation=scalars(use_session).orientation;
centroid=scalars(use_session).centroid;
%%
% threshold data

frame_data=single(depth_bounded_rotated);
frame_data(log(depth_bounded_cable_mask_rotated)<-13)=nan;
frame_data(depth_bounded_rotated<10)=0;
frame_data=reshape(frame_data,80^2,[]);

frame_data=frame_data';
replace_idx=isnan(frame_data);

for i=1:size(frame_data,1)
  idx=isnan(frame_data(i,:));
  this_frame=frame_data(i,:);
  frame_data(i,idx)=mean(this_frame(this_frame>0));
end

iters=5;

for i=1:iters
  i
  mu=mean(frame_data,1);
  scores=bsxfun(@minus,frame_data,mu)*features.pca.coeffs(:,1:50);
  recon=scores(:,1:20)*features.pca.coeffs(:,1:20)';
  recon=bsxfun(@plus,recon,mu);
  frame_data(replace_idx)=recon(replace_idx);
end

%% projections we want to plot

rps=kinect_gaussproj(frame_data,1e3);
rps=zscore(zscore(rps)');
[rps,delta_score,changepoints,changethresh]=kinect_analysis_changepoints(frame_data,'delta_thresh',.05,'smooth_sig',2);
frame_data=reshape(frame_data',80,80,[]);
recon=reshape(recon',80,80,[]);
spines=squeeze(mean(frame_data(40:60,:,:)));

%%




%%
% cross
max_lag=15;
gcamp_dff=photometry{use_session}.kin.proc.data(:,1)./...
  photometry{use_session}.kin.proc.data_baseline(:,1);
rcamp_dff=photometry{use_session}.kin.proc.data(:,2)./...
  photometry{use_session}.kin.proc.data_baseline(:,2);

[signal1_win,t_win]=markolab_vec2mat(zscore(gcamp_dff),20,19);
signal2_win=markolab_vec2mat(zscore(rcamp_dff),20,19);

rmat=zeros(max_lag*2+1,numel(t_win));

for i=1:numel(t_win)
    rmat(:,i)=xcorr(signal1_win(:,i),signal2_win(:,i),max_lag,'coeff');
end

%%

max_lag_pc=90;
rmat_pcs_gcamp=zeros(max_lag_pc*2+1,30);
rmat_pcs_rcamp=zeros(size(rmat_pcs_gcamp));

for i=1:size(rmat_pcs_gcamp,2)
   rmat_pcs_gcamp(:,i)=xcorr(gcamp_dff,scores(:,i),max_lag_pc,'coeff');
   rmat_pcs_rcamp(:,i)=xcorr(rcamp_dff,scores(:,i),max_lag_pc,'coeff');
end

%% panel for the bob-eo

use_ts=metadata(use_session).depth_timestamps_kinect;
use_ts=use_ts-min(use_ts);
use_ts=use_ts/1e3;
smooth_tau=7;
lags=[-max_lag:max_lag]./30;

ax=[];
fig=figure();
ax(1)=subplot(4,1,1);
imagesc(use_ts,[],imgaussfilt(rps,1));
axis xy;
ylabel('Rnd proj.');

colormap(bone);
caxis([-1 1]);
freezeColors();
box off;
axis off;

ax(2)=subplot(4,1,2);
imagesc(use_ts(t_win+1),lags,rmat);
colormap(jet)
caxis([-1 1]);
ylim([-max_lag/30 max_lag/30]);
ylabel('Lag (s)');
freezeColors();
box off;
set(gca,'TickDir','out','YTick',[-max_lag/30 max_lag/30],'xtick',[]);

ax(3)=subplot(4,1,3);
plot(use_ts(t_win+1),rmat(lags==0,:),'w-');
ylim([-1 1]);
ylabel('Correlation');
box off;
set(gca,'TickDir','out','YTick',[-1 1],'xtick',[]);

gcamp_dff=(gcamp_dff-min(gcamp_dff(100:end)))./(max(gcamp_dff(100:end))-min(gcamp_dff(100:end)));
rcamp_dff=(rcamp_dff-min(rcamp_dff(100:end)))./(max(rcamp_dff(100:end))-min(rcamp_dff(100:end)));

ax(4)=subplot(4,1,4);
plot(use_ts,conv(gcamp_dff,ones(smooth_tau,1)./smooth_tau,'same'),'g-');
hold on;
plot(use_ts,conv(rcamp_dff,ones(smooth_tau,1)./smooth_tau,'same'),'r-');
box off;
ylim([0 1]);
axis off;

linkaxes(ax,'x');
xlim([430 494]);

xlims=xlim;
ylims=ylim;

h=line([xlims(1)-2 xlims(1)-2],[ylims(1) ylims(1)+1],'color','w','parent',ax(4));
h2=line([xlims(1)-2 xlims(1)+3],[ylims(1) ylims(1)],'color','w','parent',ax(4));

set(h,'clipping','off');
set(h2,'clipping','off');
set(fig,'InvertHardCopy','off');

%% panels related to xcorr

[val,idx]=max(rmat_pcs_gcamp);
[~,idx2]=sort(idx);


figure();
ax(1)=subplot(2,2,1);
imagesc([-max_lag_pc:max_lag_pc]/30,[],zscore(rmat_pcs_gcamp(:,idx2))')
colormap(jet);
caxis([-3 3]);
xlim([-max_lag_pc max_lag_pc]/30);
ylim([.5 30.5]);
set(gca,'XTick',[-max_lag_pc max_lag_pc]/30,'YTick',[0.5 30.5],'YTickLabel',[1 30]);

ax(2)=subplot(2,2,2);
imagesc([-max_lag_pc:max_lag_pc]/30,[],zscore(rmat_pcs_rcamp(:,idx2))')
xlim([-max_lag_pc max_lag_pc]/30);
ylim([.5 30.5]);
set(gca,'XTick',[-max_lag_pc max_lag_pc]/30,'YTick',[0.5 30.5],'YTickLabel',[1 30]);
colormap(jet)
caxis([-3 3]);

[val,idx]=max(rmat_pcs_rcamp);
[~,idx2]=sort(idx);

ax(3)=subplot(2,2,3);
imagesc([-max_lag_pc:max_lag_pc]/30,[],zscore(rmat_pcs_gcamp(:,idx2))')
colormap(jet);
caxis([-3 3]);
xlim([-max_lag_pc max_lag_pc]/30);
ylim([.5 30.5]);
set(gca,'XTick',[-max_lag_pc max_lag_pc]/30,'YTick',[0.5 30.5],'YTickLabel',[1 30]);

ax(4)=subplot(2,2,4);
imagesc([-max_lag_pc:max_lag_pc]/30,[],zscore(rmat_pcs_rcamp(:,idx2))')
colormap(jet)
caxis([-3 3]);
xlim([-max_lag_pc max_lag_pc]/30);
ylim([.5 30.5]);
set(gca,'XTick',[-max_lag_pc max_lag_pc]/30,'YTick',[0.5 30.5],'YTickLabel',[1 30]);

% use the reconstruciton to do stuff

%%

fig=figure();
[r,lags]=xcorr(zscore(gcamp_dff),zscore(delta_score),90,'coeff');
[r2]=xcorr(zscore(rcamp_dff),zscore(delta_score),90,'coeff');

plot(lags/30,r,'g-');
hold on
plot(lags/30,r2,'r-');



%
% f%%
% % now either put recon back in
%
% full_recon=zeros(424,512,size(frame_data,3),'int16');
% box_size=[80 80];
%
% for i=1:size(frame_data,3)
%
  insert_mouse=recon(:,:,i);
%   insert_mouse=imgaussfilt(imrotate(insert_mouse,orientation(i),'bilinear','crop'),.1);
%
%
%   coords_x=centroid(i,1)-(box_size(2)/2-1):centroid(i,1)+box_size(2)/2;
%   coords_y=centroid(i,2)-(box_size(1)/2-1):centroid(i,2)+box_size(1)/2;
%
%   full_recon(coords_y,coords_x,i)=insert_mouse;
%   figure(1);
%   imagesc(full_recon(:,:,i));caxis([10 80]);
%   pause(.01);
%
% end
%
%
% %%
