%% simple example of labeling
%
%
%
%


load('experiment_data_rps.mat','rps')
load('arhmm_labels.mat','state_labels','filenames');
load('experiment_data_scores.mat','frame_idx');
load('experiment_data_neural.mat','photometry');
load('use_session.mat','use_session');
load('phase_randomized_photometry.mat','phase_rnds');

%%

use_photometry=photometry(use_session);
use_rps=rps(use_session);
plot_session=2;
%use_labels=state_labels{use_session}(~isnan(frame_idx{use_session}));
camera_fs=30;
% 15781 session 9 10659    11309

%%

[states durations usage syllable_idx starts stops]=kinect_syllable_durations(state_labels);

% sort by usage

[~,usage_idx]=sort(usage,'descend');

%%

% plot a simple example with spines, rps, photometry, labels...

plot_labels=state_labels{use_session(plot_session)}(~isnan(frame_idx{use_session(plot_session)}));
plot_rps=zscore(zscore(use_rps{plot_session})');
plot_photometry=use_photometry{plot_session};

figs.model_example=figure('PaperPositionMode','auto','position',[100 100 600 600]);
tvec=[1:size(use_rps,2)]/camera_fs;
ax(1)=subplot(2,1,1);imagesc(tvec,[],imgaussfilt(use_rps,1));axis off;
caxis([-1 1]);colormap(bone);
freezeColors();
ax(2)=subplot(2,1,2);imagesc(tvec,[],use_labels);axis off;
linkaxes(ax,'x');
xlim([355 376]);
h=line([355 357],[1.6 1.6],'color','w','parent',ax(2));
colormap(jet);
freezeColors();
%whitebg(figs.model_example);
set(figs.model_example,'color',[0 0 0],'InvertHardCopy','off');
set(h,'clipping','off');
markolab_multi_fig_save(figs.model_example,'~/Desktop/quickfigs',[ 'systemclub_modeling_example'] ,'eps,png,fig','renderer','painters');


%%
% analyze modeling same as pca, can likely copy paste code here (onsets and
% offsets???)
% loop through the 10 most commonly occurring syllables, take xcorr
%
%
syllable_win={};
photometry_win={};
photometry_rnd_win={};
durations={};
win_size=120;
nrands=1e3;
upd=kinect_proctimer(length(use_session)*45);
counter=1;
 
for i=1:length(use_session)
    for j=1:45
        tmp=syllable_idx{usage_idx(j)}{use_session(i)}(~isnan(frame_idx{use_session(i)}));
        % get onsets
        newidx=[1:length(tmp)-1];
        onsets=[0;~tmp(newidx)&tmp(newidx+1)];
        offsets=[0;tmp(newidx)&~tmp(newidx+1)];
        dur=find(offsets)-find(onsets);
        [syllable_win{j,i},~,to_keep]=markolab_win_data(find(onsets),120,tmp);
        photometry_win{j,i}=markolab_win_data(find(onsets),120,use_photometry{i}.kin.ref.data);
        upd(counter);
        counter=counter+1;
        durations{j,i}=dur(to_keep);
    end
end

% take average peak value, use this for a first-pass
%%
ave_win=zeros(win_size*2+1,45);

for i=1:45
    tmp=cat(2,photometry_win{i,:});
    ave_win(:,i)=mean(zscore(tmp)');
end

ave_rnds=zeros(win_size*2+1,45,nrands);
photometry_rnd_win={};
upd=kinect_proctimer(nrands*45);
counter=1;

for ii=1:nrands
    
    for i=1:45
        
        for j=1:length(use_session)
            
            tmp=syllable_idx{usage_idx(i)}{use_session(j)}(~isnan(frame_idx{use_session(j)}));
            % get onsets
            newidx=[1:length(tmp)-1];
            onsets=[0;~tmp(newidx)&tmp(newidx+1)];
            photometry_rnd_win{i,j}=markolab_win_data(find(onsets),120,phase_rnds{j}.ref(:,ii));
            
            
        end
        
        tmp=cat(2,photometry_rnd_win{i,:});
        ave_rnds(:,i,ii)=mean(zscore(tmp)');
        upd(counter)
        counter=counter+1;
    end
    
end

save('syllable_id_corr_phasernd.mat','ave_rnds','ave_win')

%%
clims=[-10 10];

p_val_right=mean(repmat(ave_win,[1 1 nrands])<ave_rnds,3);
p_val_left=mean(repmat(ave_win,[1 1 nrands])>ave_rnds,3);

high_hits=find(sum(p_val_right(win_size-10:win_size+10,:)<.001)>0);
low_hits=find(sum(p_val_left(win_size-10:win_size+10,:)<.001)>0);
low_hits=setdiff(low_hits,high_hits)
mu_rnds=mean(ave_rnds,3);
std_rnds=std(ave_rnds,[],3);

% 
% scales=[-10 10];

ave_win_z=(ave_win-mu_rnds)./std_rnds;

% sort the high and low hits by something interesting...
[~,tmp]=max(ave_win_z(win_size-20:win_size+20,high_hits));
[~,high_idx]=sort(tmp,'ascend');
[~,tmp]=min(ave_win_z(win_size-20:win_size+20,low_hits));
[~,low_idx]=sort(tmp,'descend');
ave_win_z=ave_win_z(:,[high_hits(high_idx)';low_hits(low_idx)']);
[val]=mean(ave_win_z(win_size:win_size+20,:));
%[~,left_idx]=sort(val,'descend');
% 
% scales_win_z=ave_win_z(:,left_idx);
% scales_win_z=(scales_win_z-(scales(1)))./(scales(2)-scales(1));
% scales_win_z(scales_win_z<0)=0;
% scales_win_z(scales_win_z>1)=1;
% scales_win_z=uint8(scales_win_z*256);
% scales_win_rgb=ind2rgb(scales_win_z,jet(256));
% 
% [r,c]=find(min(p_val_right,p_val_left)>.01);
% 
% for i=1:length(r)
%    %ms_image_rgb(r(i),c(i),:)=.01*ms_image_rgb(r(i),c(i),:); 
%    scales_win_rgb(r(i),c(i),:)=[.8 .8 .8];
% end
% 
% figure();image(scales_win_rgb);
% pause();


figs.syllable_corr=figure('PaperPositionMode','auto','position',[100 100 400 400]);
whitebg(figs.syllable_corr);
% ax(1)=subplot(1,3,1);

tvec=[-win_size:win_size]/camera_fs;
ave_win_z(abs(ave_win_z)<1)=0;
imagesc(tvec,[],ave_win_z');
colormap(hot);
set(gca,'YTick',[],'XTick',[-4:1:4],'XLim',[-3 3],'TickDir','out','TickLength',[.015 .015]);
box off;
set(figs.syllable_corr,'color',[0 0 0]);
caxis([-10 10]);

% ax(2)=subplot(1,3,2);
% [val,idx]=max(ave_win_z(100:140,high_hits));
% [~,idx2]=sort(idx,'ascend');
% imagesc(tvec,[],ave_win_z(:,high_hits(idx2))');box off;
% set(gca,'YTick',[],'XTick',[-4:1:4],'XLim',[-3 3],'TickDir','out','TickLength',[.015 .015]);
% 
% caxis([-10 10]);
% 
% 
% ax(3)=subplot(1,3,3);
% [val,idx]=min(ave_win_z(100:140,low_hits));
% [~,idx2]=sort(idx,'ascend');
% imagesc(tvec,[],ave_win_z(:,low_hits(idx2))');box off;
% set(gca,'YTick',[],'XTick',[-4:1:4],'XLim',[-3 3],'TickDir','out','TickLength',[.015 .015]);
% caxis([-10 10]);
set(figs.syllable_corr,'InvertHardcopy','off');
markolab_multi_fig_save(figs.syllable_corr,'~/Desktop/quickfigs',[ 'systems_club_syllable_corr_significant_only_hot'] ,'eps,png,fig','renderer','painters');


%%

% two color shit

syllable_win={};
photometry_win_green={};
photometry_rnd_win={};
durations={};
win_size=120;
nrands=1e3;
upd=kinect_proctimer(length(use_session)*45);
counter=1;
 
for i=1:length(use_session)
    for j=1:10
        tmp=syllable_idx{usage_idx(j)}{use_session(i)}(~isnan(frame_idx{use_session(i)}));
        % get onsets
        newidx=[1:length(tmp)-1];
        onsets=[0;~tmp(newidx)&tmp(newidx+1)];
        offsets=[0;tmp(newidx)&~tmp(newidx+1)];
        dur=find(offsets)-find(onsets);
        [syllable_win{j,i},~,to_keep]=markolab_win_data(find(onsets),120,tmp);
        photometry_win_green{j,i}=markolab_win_data(find(onsets),120,use_photometry{i}.kin.proc.data(:,1));
        photometry_win_red{j,i}=markolab_win_data(find(onsets),120,use_photometry{i}.kin.proc.data(:,2));
        upd(counter);
        counter=counter+1;
        durations{j,i}=dur(to_keep);
    end
end


ave_win_green=zeros(win_size*2+1,10);
ave_win_red=zeros(win_size*2+1,10);


for i=1:10
    tmp=cat(2,photometry_win_green{i,:});
    ave_win_green(:,i)=mean(zscore(tmp)');
    tmp=cat(2,photometry_win_red{i,:});
    ave_win_red(:,i)=mean(zscore(tmp)');
end


% for i=1:45
%    corr_vec=syllable_idx{usage_idx(i)}{use_session}(~isnan(frame_idx{use_session}));
%    newidx=[1:length(corr_vec)-1];
%    corr_vec=[0;~corr_vec(newidx)&corr_vec(newidx+1)];
%    corr_vec=conv(corr_vec,normpdf([-10:10],0,10),'same');
%    [r(i,:),lags]=xcorr(zscore(use_photometry.kin.ref.data),zscore(corr_vec),90,'coeff');
% end