%% simple example of labeling
%
%
%
%


load('experiment_data_rps.mat','rps')
load('arhmm_labels.mat','state_labels','filenames');
load('experiment_data_scores.mat','frame_idx');
load('experiment_data_neural.mat','photometry');

%%

use_session=9;
use_photometry=photometry{use_session};
use_rps=zscore(zscore(rps{use_session})');
use_labels=state_labels{use_session}(~isnan(frame_idx{use_session}));
camera_fs=30;
% 15781 session 9 10659    11309

%%

[states durations usage syllable_idx starts stops]=kinect_syllable_durations(state_labels);

% sort by usage

[~,usage_idx]=sort(usage,'descend');

%%

% plot a simple example with spines, rps, photometry, labels...

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

for i=1:45
   corr_vec=syllable_idx{usage_idx(i)}{use_session}(~isnan(frame_idx{use_session}));
   newidx=[1:length(corr_vec)-1];
   corr_vec=[0;~corr_vec(newidx)&corr_vec(newidx+1)];
   corr_vec=conv(corr_vec,normpdf([-10:10],0,10),'same');
   [r(i,:),lags]=xcorr(zscore(use_photometry.kin.ref.data),zscore(corr_vec),90,'coeff');
end