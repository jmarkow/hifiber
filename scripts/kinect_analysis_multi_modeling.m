% then for all syllables, get significance, sort by z-score

beh=extract_object.get_behavior_model;
phot=extract_object.get_photometry;
usage=beh.get_syllable_usage;
usage=usage/sum(usage);

[~,idx]=sort(usage,'descend');
nsyllables=sum(usage>.001);

for ii=1:length(extract_object)
    
    num=ii;

    if length(phot(num).traces)==0
        continue;
    end
    
    nrands=phot(num).options.photometry.nrands;
     
    max_lag=90;
    
    syll_corr_gcamp=zeros(max_lag*2+1,nsyllables);
    syll_corr_rcamp=zeros(max_lag*2+1,nsyllables);

    syll_corr_gcamp_rnd=zeros(max_lag*2+1,nsyllables,nrands);
    syll_corr_rcamp_rnd=zeros(max_lag*2+1,nsyllables,nrands);

    kernel=normpdf(-10:10,0,1.5);
    timer_upd=kinect_extract.proc_timer(nsyllables);
    
    for i=1:nsyllables

%         onset_vec=zeros(size(extract_object(num).behavior_model.labels));
%         onset_vec(extract_object(num).behavior_model.state_starts{usage_idx(i)})=1;
%         onset_vec=conv(onset_vec,kernel,'same');

        % swap out with triggered average
        
        matches=beh(num).state_starts{idx(i)};
        %matches=extract_object(num).template_match(idx(i),1);
        
        if isempty(matches)
            continue;
        end
        
        syll_corr_gcamp(:,i)=nanmean(kinect_extract.window_data(zscore(phot(num).traces(1).dff),...
            matches,...    
            max_lag),2);
        syll_corr_rcamp(:,i)=nanmean(kinect_extract.window_data(zscore(phot(num).traces(2).dff),...
            matches,...    
            max_lag),2);

        for j=1:nrands       
            syll_corr_gcamp_rnd(:,i,j)=nanmean(kinect_extract.window_data(zscore(phot(num).traces(1).dff_rnd(:,j)),...
                matches,...
                max_lag),2);
            syll_corr_rcamp_rnd(:,i,j)=nanmean(kinect_extract.window_data(zscore(phot(num).traces(2).dff_rnd(:,j)),...
                matches,...    
                max_lag),2);
        end

        timer_upd(i);

    end
    
    extract_object(ii).neural_data.analysis.syll_corr_gcamp=syll_corr_gcamp;
    extract_object(ii).neural_data.analysis.syll_corr_rcamp=syll_corr_rcamp;
    extract_object(ii).neural_data.analysis.syll_corr_gcamp_rnd=syll_corr_gcamp_rnd;
    extract_object(ii).neural_data.analysis.syll_corr_rcamp_rnd=syll_corr_gcamp_rnd;

    
end

% do the same for pcs

