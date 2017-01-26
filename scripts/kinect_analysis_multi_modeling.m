%%%% script for multi-color modeling analysis


% like for the talk we should align both to rp transitions first

% then for all syllables, get significance, sort by z-score

num=8;
nrands=phot(1).options.photometry.nrands;
usage=extract_object.get_syllable_usage;
[~,usage_idx]=sort(usage,'descend');
max_lag=90;
nsyllables=40;
corr_mat_gcamp=zeros(max_lag*2+1,nsyllables);
corr_mat_rcamp=zeros(max_lag*2+1,nsyllables);

corr_mat_gcamp_rnd=zeros(max_lag*2+1,nsyllables,nrands);
corr_mat_rcamp_rnd=zeros(max_lag*2+1,nsyllables,nrands);

kernel=normpdf(-10:10,0,3);
timer_upd=kinect_extract.proc_timer(nsyllables);
% obj=extract_object(num);
% obj.neural_data.photometry.invert;
% obj.neural_data.photometry.get_baseline;
% obj.neural_data.photometry.subtract_baseline;

for i=1:nsyllables

    onset_vec=zeros(size(extract_object(num).behavior_model.labels));
    onset_vec(extract_object(num).behavior_model.state_starts{usage_idx(i)})=1;
    onset_vec=conv(onset_vec,kernel,'same');
    
    corr_mat_gcamp(:,i)=xcorr(zscore(onset_vec),zscore(phot(num).traces(1).dff),max_lag,'coeff');
    corr_mat_rcamp(:,i)=xcorr(zscore(onset_vec),zscore(phot(num).traces(2).dff),max_lag,'coeff');
    
    for j=1:nrands       
        corr_mat_gcamp_rnd(:,i,j)=xcorr(zscore(onset_vec),zscore(phot(num).traces(1).dff_rnd(:,j)),max_lag,'coeff');
        corr_mat_rcamp_rnd(:,i,j)=xcorr(zscore(onset_vec),zscore(phot(num).traces(2).dff_rnd(:,j)),max_lag,'coeff');
    end
    
    timer_upd(i);

end

% do the same for pcs

