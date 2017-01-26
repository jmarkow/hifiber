ax2=axes('ydir','rev','units','pixels','position',[650 474-160 160 160]);
h2=imagesc(rand(size(depth_bounded_rotated(:,:,1))),'parent',ax2);caxis([0 40]);
axis(ax2,'off');

ax3=axes('units','pixels','position',[520 200 400 100]);
h3=imagesc(disp_rps(:,frames(1):frames(1)+rps_width));
h4=patch([0 rps_width+1 rps_width+1 0],[ 1 1 size(rps,1) size(rps,1) ],0,'facecolor',[0 0 0],'edgecolor','none');
caxis([-1 1]);
axis(ax3,'off');
set(ax3,'xlim',[1 rps_width]);

score_plot=filter(ones(5,1)/5,1,zscore(scores(:,:)))+repmat([size(scores,2):-1:1]*2,[size(scores,1) 1]);

ax4=axes('units','pixels','position',[520 70 400 100]);
h5=plot(score_plot(frames(1):frames(1)+rps_width,1:10),'linewidth',1.5);
ylimits=ylim(ax4);
h6=patch([0 rps_width+1 rps_width+1 0],[ ylimits(1) ylimits(1) ylimits(2) ylimits(2) ],0,'facecolor',[0 0 0],'edgecolor','none');
axis(ax4,'off');
set(ax4,'xlim',[1 rps_width]);

%open(v);
timer_upd=kinect_extract.proc_timer(length(frames));

for i=frames
    set(h,'XData',xx,'YData',yy,'ZData',depth_masked(:,:,i));
    set(h2,'CData',depth_bounded_rotated(:,:,i));
    if (i-(frames(1)-1))<=rps_width
        set(h4,'xdata',[i-(frames(1)-1) rps_width+1 rps_width+1 i-(frames(1)-1)]');
        set(h6,'xdata',[i-(frames(1)-1) rps_width+1 rps_width+1 i-(frames(1)-1)]');
    else
        set(h3,'cdata',disp_rps(:,i-rps_width:i));
        for j=1:length(h5)
            h5(j).YData=score_plot(i-rps_width:i,j);
        end
    end
    view(ax,az(i-(frames(1)-1)),el(i-(frames(1)-1)));
    %im=getframe(fig);
    %writeVideo(v,im.cdata);
    pause(eps);
    timer_upd(i-(frames(1)-1));
end

%close(v)

%% two color z-score analysis (take corr_mat)

corr_mat_gcamp_mu=squeeze(mean(corr_mat_gcamp_rnd,3));
corr_mat_gcamp_std=squeeze(std(corr_mat_gcamp_rnd,[],3));
corr_mat_rcamp_mu=squeeze(mean(corr_mat_rcamp_rnd,3));
corr_mat_rcamp_std=squeeze(std(corr_mat_rcamp_rnd,[],3));

corr_mat_gcampz=(corr_mat_gcamp-corr_mat_gcamp_mu)./corr_mat_gcamp_std;
corr_mat_rcampz=(corr_mat_rcamp-corr_mat_rcamp_mu)./corr_mat_rcamp_std;

% p-val (right tail)

use_idx=40:140;
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

pos_hits=find(pval_right_rcamp<.05&~isnan(min(corr_mat_rcamp))');
neg_hits=find(pval_left_rcamp<.05&~isnan(min(corr_mat_rcamp))');

[val,idx]=max(corr_mat_rcampz(use_idx,pos_hits));
[~,pos_idx2]=sort(idx,'ascend');

[val,idx]=min(corr_mat_rcampz(use_idx,neg_hits));
[~,neg_idx2]=sort(idx,'descend');





