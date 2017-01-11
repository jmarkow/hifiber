%%%% script for multi-color modeling analysis


% like for the talk we should align both to rp transitions first

% then for all syllables, get significance, sort by z-score

num=4;

usage=cellfun(@length,extract_object(num).behavior_model.state_durations);
[~,idx]=sort(usage,'descend');
max_lag=90;
nsyllables=50;
corr_mat=zeros(max_lag*2+1,nsyllables);
corr_mat2=zeros(max_lag*2+1,nsyllables);
kernel=normpdf(-10:10,0,2);
timer_upd=kinect_extract.proc_timer(nsyllables);

for i=1:nsyllables

    onset_vec=zeros(size(extract_object(num).behavior_model.labels));
    onset_vec(extract_object(num).behavior_model.state_starts{idx(i)})=1;
    onset_vec=conv(onset_vec,kernel,'same');
    corr_mat(:,i)=xcorr(zscore(onset_vec),zscore(phot.traces(1).baseline_rem),max_lag,'coeff');
    corr_mat2(:,i)=xcorr(zscore(onset_vec),zscore(phot.traces(2).baseline_rem),max_lag,'coeff');
    timer_upd(i);

end

% do the same for pcs
