%
%
%
%
%
%
%

% changepoints stuff

% get the std from the bootstrap

% maybe just insert a loop for this


camera_fs=30;
conditions=fieldnames(lag_mat);
conditions(strcmp(conditions,'delta_model'))=[];
conditions(strcmp(conditions,'delta_bin_midpoints'))=[];
mouse='15784';
mu={};
ymin={};
ymax={};

for i=1:length(conditions)
   mu{i}=mean(lag_mat.(conditions{i}).ref);
   ci=bootci(1e3,{@mean,lag_mat.(conditions{i}).ref},'type','per','alpha',.01);
   ymin{i}=ci(1,:)';
   ymax{i}=ci(2,:)';
   bootstr_mu{i}=mean(rnd_summary.(conditions{i}).ref);
   bootstr_min{i}=-3*std(rnd_summary.(conditions{i}).ref)';
   bootstr_max{i}=3*std(rnd_summary.(conditions{i}).ref)';
end

figs.changepoints=figure('PaperPositionMode','auto','Position',[100 100 800 250]);
whitebg(figs.changepoints);

for i=1:length(conditions)
    g(1,i)=gramm('x',([-90:90]/camera_fs)','y',{mu{i};bootstr_mu{i}},'ymin',{ymin{i};bootstr_min{i}},'ymax',{ymax{i};bootstr_max{i}},'color',1:2);
    g(1,i).axe_property('ylim',[-.05 .25],'xlim',[-3 3],'xtick',[-3 0 3],'ytick',[-.05:.05:.25],'TickDir','out','TickLength',[.015 .015]);
    %g(1,i).geom_line;
    g(1,i).geom_interval('geom','solid_area')
    g(1,i).no_legend;
    g(1,i).set_color_options('map',[1 0 0;.8 .8 .8],'chroma',50);
    g(1,i).set_title(conditions{i})
    g(1,i).set_names('x','','y','');
end

g.draw
% 
% for i=1:length(g)
%     offsetAxes(g(i).facet_axes_handles);
% end
set(figs.changepoints,'color',[0 0 0]);
set(figs.changepoints,'color',[0 0 0],'InvertHardcopy','off')

% pdf is the only output that looks half decent here

markolab_multi_fig_save(figs.changepoints,'~/Desktop/quickfigs',[ mouse '_systemclub_changepoint_xcorr'] ,'eps,png,fig','renderer','painters');

%%

files={};
files{1}='/Users/jmarkow/Desktop/workspace/photometry/15781/analysis/changepoint_analysis.mat';
files{2}='/Users/jmarkow/Desktop/workspace/photometry/15783/analysis/changepoint_analysis.mat';
files{3}='/Users/jmarkow/Desktop/workspace/photometry/15784/analysis/changepoint_analysis.mat';


load(files{1});

conditions=fieldnames(lag_mat);
conditions(strcmp(conditions,'delta_model'))=[];
conditions(strcmp(conditions,'delta_bin_midpoints'))=[];

for i=2:length(files)
   new_data=load(files{2});
   for j=1:length(conditions)
       lag_mat.(conditions{j}).ref=[lag_mat.(conditions{j}).ref;new_data.lag_mat.(conditions{j}).ref];
       rnd_summary.(conditions{j}).ref=[rnd_summary.(conditions{j}).ref;new_data.rnd_summary.(conditions{j}).ref];
   end
end
figs.d1_summary=figure('PaperPositionMode','auto','position',[100 100 500 200]);
whitebg(figs.d1_summary);

for i=1:length(conditions)
   rnd_max=max(rnd_summary.(conditions{i}).ref');
   
   lag_max=max(lag_mat.(conditions{i}).ref');
   rnd_max=randsample(rnd_max,length(lag_max),false);
   category=[ones(1,length(rnd_max)) ones(1,length(lag_max))*2];
   g(1,i)=gramm('x',[ones(1,length(rnd_max)) ones(1,length(lag_max))*2],'y',[rnd_max(:)' lag_max(:)'],'color',category);
   g(1,i).geom_jitter();
   g(1,i).no_legend();
   g(1,i).set_names('x','','y','');
   g(1,i).axe_property('ylim',[-.2 .7],'TickDir','out','TickLength',[.015 .015],'YTick',[-.2:.2:.8]); 
end

g.draw;
set(figs.d1_summary,'color',[0 0 0],'InvertHardcopy','off');
markolab_multi_fig_save(figs.d1_summary,'~/Desktop/quickfigs',[ 'systemsclub_d1summary'] ,'eps,png,fig','renderer','painters');

%%

% combine across mice?

files={};
files{1}='/Users/jmarkow/Desktop/workspace/photometry/15781/analysis/timescale_analysis.mat';
files{2}='/Users/jmarkow/Desktop/workspace/photometry/15783/analysis/timescale_analysis.mat';
files{3}='/Users/jmarkow/Desktop/workspace/photometry/15784/analysis/timescale_analysis.mat';

max_r_all={};
max_r_rnd={};
tmp1=[];

for i=1:length(files)
    load(files{i},'max_r','score_pks','rnd_max_r');
    max_r_all=[max_r_all max_r];
    max_r_rnd=[max_r_rnd rnd_max_r];
    tmp1=[tmp1;cellfun(@(x) median(diff(x)),score_pks)];
end

max_r_rnd=cat(1,max_r_rnd{:});
max_r_all=cat(1,max_r_all{:});
mins=min(max_r_all,[],2);
maxs=max(max_r_all,[],2);

max_r_all_raw=(max_r_all-repmat(mins,[1 size(max_r_all,2)]))./(repmat(maxs,[1 size(max_r_all,2)])-repmat(mins,[1 size(max_r_all,2)]))
max_r_all=max_r_all./squeeze(std(max_r_rnd,[],3));

figs.timescale=figure('PaperPositionMode','auto','Position',[100 100 600 250]);
whitebg(figs.timescale);

tmp1=median(tmp1)*1e3;
tmp2=mean(max_r_all);
tmp3=mean(max_r_all_raw);
newx=1:.1:length(tmp1);
tmp1=interp1(1:length(tmp1),tmp1,newx,'spline');
tmp2=interp1(1:length(tmp2),tmp2,newx,'spline');
tmp3=interp1(1:length(tmp3),tmp3,newx,'spline');
ci=bootci(1e3,{@mean,max_r_all},'type','per','alpha',.01);
ci2=bootci(1e3,{@mean,max_r_all_raw},'type','per','alpha',.01);
ci=interp1(1:size(ci,2),ci',1:.1:size(ci,2),'spline')';
ci2=interp1(1:size(ci2,2),ci2',1:.1:size(ci2,2),'spline')';

g(1,1)=gramm('x',tmp1*1/30,'y',tmp3,'ymin',ci2(1,:),'ymax',ci2(2,:));
g(1,1).geom_interval('geom','solid_area');
g(1,1).no_legend;
g(1,1).set_names('x','','y','');
g(1,1).axe_property('ylim',[0 1],'ytick',[0:.5:1],'TickDir','out','TickLength',[.015 .015],'xlim',[0 1400],'xtick',[0:200:1400])
g(1,2)=gramm('x',tmp1*1/30,'y',tmp2,'ymin',ci(1,:),'ymax',ci(2,:));
%g.geom_line;
g(1,2).geom_interval('geom','solid_area');
g(1,2).no_legend;
g(1,2).set_names('x','','y','');
g(1,2).axe_property('ylim',[1.5 3.5],'ytick',[0:.5:4],'TickDir','out','TickLength',[.015 .015],'xlim',[0 1400],'xtick',[0:200:1400])


g.draw
%offsetAxes(g.facet_axes_handles);

set(figs.timescale,'color',[0 0 0],'InvertHardcopy','off');
markolab_multi_fig_save(figs.timescale,'~/Desktop/quickfigs',['systemsclub_timescale'] ,'eps,png,fig','renderer','painters');

%%

% random projection example, use spines too??? from 15783

load('experiment_data_neural.mat','photometry');
load('experiment_data_metadata.mat','metadata');
load('experiment_data_rps.mat','rps');
load('experiment_data_pca.mat','features');
load('experiment_data_scores.mat','scores','frame_idx');

use_session=6;
use_photometry=photometry{use_session};
use_metadata=metadata(use_session);
%use_delta_score=delta_score{use_session};
use_rps=zscore(zscore(rps{use_session})');
%use_scores=scores{use_session}(:,~isnan(frame_idx{use_session}));
load(use_metadata.frames_file,'depth_bounded_rotated','depth_bounded_cable_mask_rotated');


delta_thresh=.1;
delta_win=3;
smooth_sig=1;

deltas=markolab_deltacoef(use_rps,3); % delta coefficients, lag of 4 frames
delta_score=sum(abs(deltas)>delta_thresh); % binarize deltas
h=normpdf([-10:10],0,smooth_sig); % gauss smooth
delta_vec=conv(delta_score,h,'same'); % convolve

%%


recon=depth_bounded_rotated.*int16(log(depth_bounded_cable_mask_rotated)>-14);
tvec=[1:size(use_rps,2)]./30;
ax=[];
figs.example_raster=figure('PaperPositionMode','Auto','position',[100 100 700 600]);
whitebg(figs.example_raster);
ax(1)=subplot(4,1,1);
imagesc(tvec,[],squeeze(mean(recon(35:45,:,:))));colormap(bone);
axis off;
axis xy;
caxis([0 60]);
ax(2)=subplot(4,1,2);
imagesc(tvec,[],imgaussfilt(use_rps,1));colormap(bone);
axis off;
axis xy;
caxis([-1 1]);
ax(3)=subplot(4,1,3);
plot(tvec,delta_vec,'w-');
set(ax(3),'ydir','rev');
ylim([400 800]);
axis off;
ax(4)=subplot(4,1,4);
plot(tvec,use_photometry.kin.ref.data,'g-');ylim([0 3]);
linkaxes(ax,'x');
xlim([220 240]);
axis off;
h=line([218 218],[-1 0],'color','w');
h2=line([218 220],[-1 -1],'color','w');
set(h,'clipping','off');
set(h2,'clipping','off');
set(figs.example_raster,'color',[0 0 0],'InvertHardCopy','off');

markolab_multi_fig_save(figs.example_raster,'~/Desktop/quickfigs',['systemsclub_rasterexample'] ,'eps,png,fig','renderer','painters');

%%

load('experiment_data_pca.mat','features');
load('pca_analysis.mat','lag_mat','rnd_summary');

zmu=zeros(length(lag_mat),size(lag_mat(1).ref,2));
pmu_right=zeros(size(zmu));
pmu_left=zeros(size(zmu));
mu=zeros(size(zmu));

for i=1:length(lag_mat)
   
   % express as X fold over the phase rnd
   
   mu(i,:)=mean(lag_mat(i).ref);
   zmu(i,:)=(mu(i,:)-mean(rnd_summary(i).ref))./std(rnd_summary(i).ref);
   
   % cut off by p-value?
   
   % exclude all significant pcs, two-tailed .05
   
   pmu_left(i,:)=mean(repmat(mu(i,:),[size(rnd_summary(i).ref,1) 1])>rnd_summary(i).ref);
   pmu_right(i,:)=mean(repmat(mu(i,:),[size(rnd_summary(i).ref,1) 1])<rnd_summary(i).ref);
      
end

% only include pcs with p-values<.01 in reasonable lags (-1/+1sec)

alpha_cut=.05;
right_idx=find(sum(pmu_right(:,30:150)'<=alpha_cut)>0);
left_idx=find(sum(pmu_left(:,30:150)'<=alpha_cut)>0);

zmu_right=zmu(right_idx,:);
zmu_left=zmu(left_idx,:);

% first make pseudocolor plots

%%

% SCHEME (average on top)
% sorted by lag, then show PC next to it (very small)

figs.pc_pseudocolor=figure('PaperPositionMode','auto','Position',[100 100 200 800]);
whitebg(figs.pc_pseudocolor);

lags=(size(zmu,2)-1)/2;
lagvec=[-lags:lags]/30;
[maxs,idx]=max(zmu_right(:,30:150)');
mins=min(zmu_right');
[~,idx2]=sort(idx,'ascend');
nsamples=size(zmu_right,2);
whitebg(figs.pc_pseudocolor);
zmu_right_scaled=(zmu_right-repmat(mins(:),[1 nsamples]))./(repmat(maxs(:),[1 nsamples])-repmat(mins(:),[1 nsamples]));

 ax=[];
ax(1)=subplot(5,1,1);
plot(lagvec,mean(zmu));box off
set(gca,'xtick',[],'ytick',[-1 2],'ylim',[-1 2],'TickDir','out','TickLength',[.015 .015],'FontSize',12);
ax(2)=subplot(5,1,2:5);
imagesc(lagvec,[],zmu_right_scaled(idx2,:));box off;

colormap(jet);
set(gca,'TickDir','out','TickLength',[.015 .015],'FontSize',12,'YTick',[]);
caxis([0 1.2]);

% 
% [maxs]=max(zmu_left');
% [mins,idx]=min(zmu_left');
% [~,idx2]=sort(idx,'ascend');
% nsamples=size(zmu_left,2);
% whitebg(figs.pc_pseudocolor);
% zmu_left_scaled=(zmu_left-repmat(mins(:),[1 nsamples]))./(repmat(maxs(:),[1 nsamples])-repmat(mins(:),[1 nsamples]));
% ax(2)=subplot(1,2,2);
% imagesc(lagvec,[],zmu_left_scaled(idx2,:));box off;
% set(gca,'TickDir','out','TickLength',[.015 .015],'FontSize',12,'YTick',[]);
% xlim([-1 1]);
 linkaxes(ax,'x');
xlim([-2 2]);
set(figs.pc_pseudocolor,'color',[0 0 0],'InvertHardcopy','off');
markolab_multi_fig_save(figs.pc_pseudocolor,'~/Desktop/quickfigs',['systemsclub_pcpseudocolor'] ,'eps,png,fig','renderer','painters');


%%

figs.pc_sortedeigenmice=figure('PaperPositionMode','auto','Position',[100 100 300 800]);
kinect_eigenmontage(features.pca.coeffs(:,right_idx(idx2)),'k',length(right_idx),'ncols',1);
set(figs.pc_sortedeigenmice,'color',[0 0 0],'InvertHardcopy','off');
% make a mousey based on d1-corr
markolab_multi_fig_save(figs.pc_sortedeigenmice,'~/Desktop/quickfigs',['systemsclub_eigenmice_pseudocolorsorting'] ,'eps,png,fig','renderer','painters');


%%
figs.pc_reconmouse=figure('PaperPositionMode','auto','Position',[100 100 300 300]);
% recon=features.pca.coeffs(:,right_idx)*max(zmu_right')';
% set(figs.pc_reconmouse,'color',[0 0 0],'InvertHardcopy','off');
% imagesc(abs(reshape(recon,80,80)))
% caxis([.25 1]);colormap(jet);

recon=abs(reshape(recon,80,80));
ms_image=(recon-0)./(.75-0);
ms_image(ms_image<0)=0;
ms_image(ms_image>1)=1;
ms_image=uint8(ms_image*256);
ms_image_rgb=ind2rgb(ms_image,jet(256));

[r,c]=find(recon<.5);

for i=1:length(r)
   %ms_image_rgb(r(i),c(i),:)=.01*ms_image_rgb(r(i),c(i),:); 
   ms_image_rgb(r(i),c(i),:)=[.8 .8 .8];
end

image(imgaussfilt(ms_image_rgb,.5));axis off;
set(figs.pc_reconmouse,'color',[0 0 0],'InvertHardcopy','off');
markolab_multi_fig_save(figs.pc_reconmouse,'~/Desktop/quickfigs',['systemsclub_pcreconmouse'] ,'eps,png,fig','renderer','painters');
