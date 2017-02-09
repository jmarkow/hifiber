

% like for the talk we should align both to rp transitions first

% then for all syllables, get significance, sort by z-score


for ii=1:length(extract_object)
    
    num=ii;

    if length(phot(num).traces)==0
        continue;
    end
    
    nrands=phot(num).options.photometry.nrands;
    usage=extract_object.get_syllable_usage;
    usage=usage/sum(usage);
    [~,usage_idx]=sort(usage,'descend');
    max_lag=90;
    nsyllables=sum(usage>.001);
    syll_corr_gcamp=zeros(max_lag*2+1,nsyllables);
    syll_corr_rcamp=zeros(max_lag*2+1,nsyllables);

    syll_corr_gcamp_rnd=zeros(max_lag*2+1,nsyllables,nrands);
    syll_corr_rcamp_rnd=zeros(max_lag*2+1,nsyllables,nrands);

    kernel=normpdf(-10:10,0,1.5);
    timer_upd=kinect_extract.proc_timer(nsyllables);
    
    for i=1:nsyllables

        onset_vec=zeros(size(extract_object(num).behavior_model.labels));
        onset_vec(extract_object(num).behavior_model.state_starts{usage_idx(i)})=1;
        onset_vec=conv(onset_vec,kernel,'same');

        % swap out with triggered average
        
        syll_corr_gcamp(:,i)=xcorr(zscore(onset_vec),zscore(phot(num).traces(1).dff),max_lag,'coeff');
        syll_corr_rcamp(:,i)=xcorr(zscore(onset_vec),zscore(phot(num).traces(2).dff),max_lag,'coeff');

        for j=1:nrands       
            syll_corr_gcamp_rnd(:,i,j)=xcorr(zscore(onset_vec),zscore(phot(num).traces(1).dff_rnd(:,j)),max_lag,'coeff');
            syll_corr_rcamp_rnd(:,i,j)=xcorr(zscore(onset_vec),zscore(phot(num).traces(2).dff_rnd(:,j)),max_lag,'coeff');
        end

        timer_upd(i);

    end
    
    extract_object(ii).neural_data.analysis.syll_corr_gcamp=syll_corr_gcamp;
    extract_object(ii).neural_data.analysis.syll_corr_rcamp=syll_corr_rcamp;
    extract_object(ii).neural_data.analysis.syll_corr_gcamp_rnd=syll_corr_gcamp_rnd;
    extract_object(ii).neural_data.analysis.syll_corr_rcamp_rnd=syll_corr_gcamp_rnd;

    
end

% do the same for pcs

