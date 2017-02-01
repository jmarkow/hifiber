% get frames 

num=4;

depth_bounded_rotated=extract_object(num).load_oriented_frames(true);

%%
rps=extract_object(num).get_original_timebase(extract_object(num).projections.rp)';
score_plot=extract_object(num).get_original_timebase(zscore(extract_object(num).projections.rp_changepoint_score));
%frame_width=diff(frame_range);
% separate movie
frame_idx=[4.2e4 4.5e4];
frames=frame_idx(1):frame_idx(2);

frame_width=500;
kernel=normpdf([-30:30],0,3);
traces(:,1)=zscore(conv(phot(num).traces(1).dff,kernel,'same'));
traces(:,2)=zscore(conv(phot(num).traces(2).dff,kernel,'same'));

%%

% static figure to find compelling example

figure();
ax(1)=subplot(4,1,1:2);imagesc(rps);
ax(2)=subplot(4,1,3);plot(zscore(score_plot));
ax(3)=subplot(4,1,4);plot(traces);

linkaxes(ax,'x');



%%

% figure();
% ax_mouse=axes('ydir','rev','units','pixels','position',[650 474-160 160 160]);
% h_mouse=imagesc(rand(size(frames(:,:,1))),'parent',ax_mouse);caxis([0 40]);
% axis(ax_mouse,'off');



fig=figure('color',[0 0 0],'InvertHardCopy','off','Position',[200 200 800 550]);

traces(frames,:)=zscore(traces(frames,:));

ax_mouse=axes('ydir','rev','units','pixels','position',[25 200 160 160]);
h_mouse=imagesc(depth_bounded_rotated(:,:,1),'parent',ax_mouse);caxis([0 40]);
colormap(bone);
axis(ax_mouse,'off');

ax1=axes('units','pixels','position',[300 350 450 150]);
h1_im=imagesc(rps(:,frame_idx(1):frame_idx(1)+frame_width),'parent',ax1);
h1_shroud=patch([0 frame_width+1 frame_width 0],[ 1 1 size(rps,1) size(rps,1) ],0,'facecolor',[0 0 0],'edgecolor','none');
caxis([-1 1]);
colormap(bone)
axis(ax1,'off');
set(ax1,'xlim',[1 frame_width]);

%score_plot=filter(ones(5,1)/5,1,zscore(scores(:,:)))+repmat([size(scores,2):-1:1]*2,[size(scores,1) 1]);

ax2=axes('units','pixels','position',[300 225 450 100]);
h2_score=plot(score_plot(frames(1):frames(1)+frame_width),'linewidth',1.5,'color','w');
set(ax2,'YDir','rev');
ylim([0 4]);
ylimits=ylim(ax2);
h2_shroud=patch([0 frame_width+1 frame_width+1 0],[ ylimits(1) ylimits(1) ylimits(2) ylimits(2) ],0,'facecolor',[0 0 0],'edgecolor','none');
axis(ax2,'off');
set(ax2,'xlim',[1 frame_width]);

ax3=axes('units','pixels','position',[300 75 450 210]);
h3_traces(1)=plot(traces(frames(1):frames(1)+frame_width,1),'linewidth',1.5,'color','g');
hold on;
h3_traces(2)=plot(traces(frames(1):frames(1)+frame_width,2),'linewidth',1.5,'color','r');
ylim([-4 7]);
ylimits=ylim(ax3);
h3_shroud=patch([0 frame_width+1 frame_width+1 0],[ ylimits(1) ylimits(1) ylimits(2) ylimits(2) ],0,'facecolor',[0 0 0],'edgecolor','none');
axis(ax3,'off');
set(ax3,'xlim',[1 frame_width]);

v=VideoWriter('twocolor_example.mp4','mpeg-4');
v.FrameRate=30;
v.Quality=100;
open(v);

timer_upd=kinect_extract.proc_timer(length(frames));

for i=frames
    %set(h,'XData',xx,'YData',yy,'ZData',depth_masked(:,:,i));
    set(h_mouse,'CData',depth_bounded_rotated(:,:,i));
    if (i-(frames(1)-1))<=frame_width
        set(h1_shroud,'xdata',[i-(frames(1)-1) frame_width+1 frame_width+1 i-(frames(1)-1)]');
        set(h2_shroud,'xdata',[i-(frames(1)-1) frame_width+1 frame_width+1 i-(frames(1)-1)]');
        set(h3_shroud,'xdata',[i-(frames(1)-1) frame_width+1 frame_width+1 i-(frames(1)-1)]');
    else
        set(h1_im,'cdata',rps(:,i-frame_width:i));
        set(h2_score,'ydata',score_plot(i-frame_width:i));
        set(h3_traces(1),'ydata',traces(i-frame_width:i,1));
        set(h3_traces(2),'ydata',traces(i-frame_width:i,2));

%         for j=1:length(h5)
%             h5(j).YData=score_plot(i-rps_width:i,j);
%         end
     end
    %view(ax,az(i-(frames(1)-1)),el(i-(frames(1)-1)));
    im=getframe(fig);
    writeVideo(v,im.cdata);
     pause(eps);
    timer_upd(i-(frames(1)-1));
end

close(v)

%% two color z-score analysis (take corr_mat)

corr_mat_gcamp_mu=squeeze(mean(corr_mat_gcamp_rnd,3));
corr_mat_gcamp_std=squeeze(std(corr_mat_gcamp_rnd,[],3));
corr_mat_rcamp_mu=squeeze(mean(corr_mat_rcamp_rnd,3));
corr_mat_rcamp_std=squeeze(std(corr_mat_rcamp_rnd,[],3));

corr_mat_gcampz=(corr_mat_gcamp-corr_mat_gcamp_mu)./corr_mat_gcamp_std;
corr_mat_rcampz=(corr_mat_rcamp-corr_mat_rcamp_mu)./corr_mat_rcamp_std;

% p-val (right tail)

use_idx=1:181;
tmp=squeeze(max(corr_mat_gcamp_rnd(use_idx,:,:)));
tmp2=repmat(max(corr_mat_gcamp(use_idx,:))',[1 1e3]);

pval_right_gcamp=mean(tmp>tmp2,2);

tmp=squeeze(max(corr_mat_rcamp_rnd(use_idx,:,:)));
tmp2=repmat(max(corr_mat_rcamp(use_idx,:))',[1 1e3]);

pval_right_rcamp=mean(tmp>tmp2,2);

% p-val (left tail)

tmp=squeeze(min(corr_mat_gcamp_rnd(use_idx,:,:)));
tmp2=repmat(min(corr_mat_gcamp(use_idx,:))',[1 1e3]);

pval_left_gcamp=mean(tmp<tmp2,2);

tmp=squeeze(min(corr_mat_rcamp_rnd(use_idx,:,:)));
tmp2=repmat(min(corr_mat_rcamp(use_idx,:))',[1 1e3]);

pval_left_rcamp=mean(tmp<tmp2,2);

%%

% plotting, make matrix like we did for gcamp, left to right pos then right
% to left neg

pos_hits_rcamp=find(pval_right_rcamp<.05&~isnan(min(corr_mat_rcamp))');
neg_hits_rcamp=find(pval_left_rcamp<.05&~isnan(min(corr_mat_rcamp))');

pos_hits_gcamp=find(pval_right_gcamp<.05&~isnan(min(corr_mat_gcamp))');
neg_hits_gcamp=find(pval_left_gcamp<.05&~isnan(min(corr_mat_gcamp))');

pos_hits_all=unique([pos_hits_rcamp;pos_hits_gcamp]);
neg_hits_all=unique([neg_hits_rcamp;neg_hits_gcamp]);

kernel=normpdf(-50:50,0,8);
use_mat_rcampz=corr_mat_rcampz(:,pos_hits_all);
use_mat_gcampz=corr_mat_gcampz(:,pos_hits_all);


for i=1:length(pos_hits_all)
    use_mat_gcampz(:,i)=conv(use_mat_gcampz(:,i),kernel,'same');
    use_mat_rcampz(:,i)=conv(use_mat_rcampz(:,i),kernel,'same');
end

%use_idx=1:181;
[val,idx]=max(use_mat_rcampz(use_idx,:));
[~,pos_rcamp_idx2]=sort(idx,'ascend');

[val,idx]=min(use_mat_rcampz(use_idx,:));
[~,neg_rcamp_idx2]=sort(idx,'descend');

[val,idx]=max(use_mat_gcampz(use_idx,:));
[~,pos_gcamp_idx2]=sort(idx,'ascend');

[val,idx]=min(use_mat_gcampz(use_idx,:));
[~,neg_gcamp_idx2]=sort(idx,'descend');

% apply sorting to each matrix

% smoothing


% and the plotting...lay out syllablen idx to match up movies
max_lag=90;
xvec=[-90:90]/30;
clims=[-2.5 2.5];

fig=figure();

ax(1)=subplot(2,2,1);
imagesc(xvec(use_idx),[],zscore(use_mat_rcampz(use_idx,pos_rcamp_idx2))');
caxis([clims]);
axis off;

ax(2)=subplot(2,2,2);
imagesc(xvec(use_idx),[],zscore(use_mat_gcampz(use_idx,pos_rcamp_idx2))');
caxis([clims]);
axis off;

ax(3)=subplot(2,2,3);
imagesc(xvec(use_idx),[],zscore(use_mat_rcampz(use_idx,pos_gcamp_idx2))');
caxis([clims]);
set(gca,'YTick',[],'XTick',[-3:1:3],'FontSize',14);
box off;

ax(4)=subplot(2,2,4);
imagesc(xvec(use_idx),[],zscore(use_mat_gcampz(use_idx,pos_gcamp_idx2))');
caxis([clims]);
set(gca,'YTick',[],'XTick',[-3:1:3],'FontSize',14);
box off;

disp([[1:length(pos_hits_all)]' pos_hits_all(pos_rcamp_idx2(:)) pos_hits_all(pos_gcamp_idx2(:))]);
whitebg(fig);
set(fig,'Color',[0 0 0],'InvertHardCopy','off');

colormap(jet);