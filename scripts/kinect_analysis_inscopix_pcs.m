
use_obj=extract_object(2);
% get rois

nrois=size(neuron_results.C,1);
ninsc_frames=size(neuron_results.C,2);

% get pcs in inscopix timebase

use_obj.load_inscopix_timestamps;
pcs=use_obj.get_original_timebase(use_obj.projections.pca);
pcs=use_obj.get_inscopix_timebase(pcs,1:ninsc_frames);


%%

npcs=size(pcs,2);
max_lag=60;

corr_mat=zeros(max_lag*2+1,npcs,nrois);

for i=1:npcs
    for j=1:nrois
        nans=isnan(neuron_results.C(j,:));
        neuron_results.C(j,nans)=0;
        corr_mat(:,i,j)=xcorr(zscore(pcs(1:30e3,i)),zscore(neuron_results.C(j,1:30e3)),max_lag,'coeff');
    end
end

%%

% cluster based on the correlation matrix

ext
changemat=squeeze(mean(win_population(80:120,:,:)));
[coeff score]=pca(zscore(changemat));

% zscore,pca, kmeans (20 looks good!)
%%

    % do this the old-fashioned way
    
clust_choice=[5:20];
bic=[];
upd=kinect_extract.proc_timer(length(clust_choice));
options=statset('MaxIter',1e3);

for i=1:length(clust_choice);
    guess=kmeans(zscore(score(:,1:10)),clust_choice(i),'replicates',5);
    model_obj{i}=fitgmdist(zscore(score(:,1:10)),clust_choice(i),'start',guess,'regularization',1e-5,'options',options);
    bic(i)=model_obj{i}.BIC;
    upd(i);
end
    
 