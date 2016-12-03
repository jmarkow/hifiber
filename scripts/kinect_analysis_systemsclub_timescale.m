
%% timescale scan, let's do this again

% need to:
% (1) reconstruct frames
% (2) get rps
% (3) get delta score
% (4) get delta score for larger windows
% (5) check xcorr across time-scales

% 15781

mouse_name='15783';
%use_session=setdiff(8:22,[16 18:19]);
use_session=setdiff(6:13,11);
use_photometry=photometry(use_session);
use_scalars=scalars(use_session);
use_metadata=metadata(use_session);
use_rps=rps(use_session);

%


%%
% delta scores

delta_win=[2:14];
smooth_sig=.25;
delta_thresh=.1;

timer_upd=kinect_proctimer(length(delta_win)*length(use_photometry));
delta_score=cell(1,length(use_photometry));

counter=1;
for i=1:length(use_photometry)

    delta_score{i}=zeros(size(use_rps{i},1),length(delta_win));
    tmp_rps=zscore(zscore(use_rps{i})');

    for j=1:length(delta_win)
        deltas=markolab_deltacoef(tmp_rps,delta_win(j)); % delta coefficients, lag of 4 frames
        delta_score{i}(:,j)=sum(abs(deltas)>delta_thresh); % binarize deltas
        h=normpdf([-10:10],0,smooth_sig); % gauss smooth
        delta_score{i}(:,j)=conv(delta_score{i}(:,j),h,'same'); % convolve
        timer_upd(counter);
        counter=counter+1;
    end
end

%%

score_pks=cell(length(use_photometry),length(delta_win));

for i=1:length(use_photometry)
    for j=1:length(delta_win)
        threshold=mean(delta_score{i}(:,j))+.5*std(delta_score{i}(:,j));
        [~,score_pks{i,j}]=findpeaks(delta_score{i}(:,j),'minpeakheight',threshold);
    end
end

%%

win_size=3600;

max_r=cell(1,length(use_photometry));
for i=1:length(use_photometry)
    mat_photometry=markolab_vec2mat(zscore(use_photometry{i}.kin.ref.data),win_size,0);
    nwins=size(mat_photometry,2);
    max_r{i}=zeros(nwins,length(delta_win));
    for j=1:length(delta_win)
        mat_score=markolab_vec2mat(zscore(delta_score{i}(:,j)),win_size,0);
        for k=1:nwins
           r=xcorr(mat_photometry(:,k),mat_score(:,k),100,'coeff'); 
           max_r{i}(k,j)=max(r);
        end
    end
end

%% 

% for plotting interpolate data, use gramm for the rest

save([ 'analysis/' mouse_name '_timescale_analysis.mat'],'max_r','score_pks','delta_score','use_session','-v7.3');

