%%%
%%
%
%

% use two metrics, (1) average correlation between rps and model
% changepoints (2) difference in distributions

[nkappas,nsessions]=size(labels);
peakthresh=.1;
peak_dist=cell(size(labels));

timer_upd=kinect_extract.proc_timer(nkappas*nsessions);
counter=0;
test=zeros(nkappas,nsessions);

for i=1:nkappas
    
    for j=1:nsessions
       
%        [~,locs]=findpeaks(extract_object(j).projections.rp_changepoint_score,'minpeakheight',peakthresh);
%        
       if ~isnan(labels{i,j})
           model_changepoints=abs(diff([-1;labels{i,j}]))>0;
       else
           continue;
       end
%        
%        peak_dist{i,j}=zeros(size(model_changepoints));
%        
%        for k=1:length(model_changepoints)
%           peak_dist{i,j}(k)=min(abs(model_changepoints(k)-locs)); 
%        end
       
       tmp=corrcoef(zscore(model_changepoints),zscore(extract_object(j).projections.rp_changepoint_score));
       test(i,j)=tmp(2,1);
        counter=counter+1;
        timer_upd(counter);
        
    end
end