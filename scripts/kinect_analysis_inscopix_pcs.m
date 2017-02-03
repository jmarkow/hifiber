
use_obj=extract_object(1);
% get rois

nrois=size(neuron_results.C,1);
ninsc_frames=size(neuron_results.C,2);

% get pcs in inscopix timebase

use_obj.load_inscopix_timestamps;
pcs=use_obj.get_original_timebase(use_obj.projections.pca);
pcs=use_obj.get_inscopix_timebase(pcs,1:ninsc_frames);
change_score=use_obj.get_original_timebase(use_obj.projections.rp_changepoint_score);
change_score=use_obj.get_inscopix_timebase(change_score,1:ninsc_frames);

%%

npcs=30;
max_lag=30;

corr_mat_pcs=zeros(max_lag*2+1,npcs,nrois);
upd=kinect_extract.proc_timer(npcs*nrois,'frequency',20);
counter=0;
use_ts=1:ninsc_frames;

for i=1:npcs
    for j=1:nrois
        nans=isnan(neuron_results.C(j,:));
        neuron_results.C(j,nans)=0;
        corr_mat_pcs(:,i,j)=xcorr(zscore(pcs(use_ts,i)),zscore(neuron_results.C(j,use_ts)),max_lag,'coeff');
        counter=counter+1;
        upd(counter);
    end
end

%%

corr_mat_rps=zeros(max_lag*2+1,nrois);
upd=kinect_extract.proc_timer(nrois,'frequency',20);
counter=0;
use_ts=1:ninsc_frames;


for j=1:nrois
    nans=isnan(neuron_results.C(j,:));
    neuron_results.C(j,nans)=0;
    corr_mat_rps(:,j)=xcorr(zscore(change_score),zscore(neuron_results.C(j,use_ts)),max_lag,'coeff');
    counter=counter+1;
    upd(counter);
end
    

%%

% cluster based on the correlation matrix

use_idx=1:61;
trim_mat=corr_mat_pcs(use_idx,:,:);
[val,idx]=max(abs(trim_mat));
val=squeeze(val);
idx=squeeze(idx);

for i=1:size(idx,1)
    for j=1:size(idx,2)
        val(i,j)=val(i,j)*sign(trim_mat(idx(i,j),i,j));
    end 
end

%clust_mat_pcs=val;
clust_mat_pcs_raw=reshape(corr_mat_pcs,(max_lag*2+1)*npcs,nrois);
[coeff_pcs clust_mat_pcs]=pca(zscore(clust_mat_pcs_raw)');
[coeff_rps clust_mat_rps]=pca(zscore(corr_mat_rps)');
% zscore,pca, kmeans (20 looks good!)
%%

    % do this the old-fashioned way
    
clust_choice=[1:30];
bic_pcs=[];
model_obj_pcs={};

upd=kinect_extract.proc_timer(length(clust_choice));
options=statset('MaxIter',1e4);

for i=1:length(clust_choice);

    guess=kmeans(clust_mat_pcs(:,1:30),clust_choice(i),'replicates',50);
    model_obj_pcs{i}=fitgmdist(clust_mat_pcs(:,1:30),clust_choice(i),'start',guess,'regularization',1e-5,'options',options);
    bic_pcs(i)=model_obj_pcs{i}.BIC;
    upd(i);
end



%%
% do this the old-fashioned way
    
clust_choice=[1:30];
bic_rps=[];
model_obj_rps={};

upd=kinect_extract.proc_timer(length(clust_choice));
options=statset('MaxIter',1e4);

for i=1:length(clust_choice)
    guess=kmeans(clust_mat_rps(:,1:30),clust_choice(i),'replicates',50);
    model_obj_rps{i}=fitgmdist(clust_mat_rps(:,1:30),clust_choice(i),'start',guess,'regularization',1e-5,'options',options);
    bic_rps(i)=model_obj_rps{i}.BIC;
    upd(i);
end

%%

% get centroids

centroids=zeros(nrois,2);
[xx,yy]=meshgrid(1:neuron_results.options.d2,1:neuron_results.options.d1);
for i=1:nrois
   mat=neuron_results.A(:,i)./sum(neuron_results.A(:,i));
   centroids(i,1)=sum(xx(:).*mat(:));
   centroids(i,2)=sum(yy(:).*mat(:)); 
end
 