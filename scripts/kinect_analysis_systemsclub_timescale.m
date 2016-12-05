
%% timescale scan, let's do this again

% need to:
% (1) reconstruct frames
% (2) get rps
% (3) get delta score
% (4) get delta score for larger windows
% (5) check xcorr across time-scales

% 15781

delta_win=[1:15];
smooth_sig=1;
delta_thresh=.1;
win_size=3600;
nrands=1e3;

load('experiment_data_neural.mat','photometry');
load('experiment_data_rps.mat','rps');
load('use_session.mat','use_session');
load('phase_randomized_photometry.mat','phase_rnds');

use_photometry=photometry(use_session);

use_rps=rps(use_session);

%

%%
% delta scores


if exist('timescale_analysis.mat','file')~=2
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
else
    fprintf('Loading delta scores.\n');
    load('timescale_analysis.mat','delta_score');
end

%%

score_pks=cell(length(use_photometry),length(delta_win));

for i=1:length(use_photometry)
    for j=1:length(delta_win)
        threshold=mean(delta_score{i}(:,j))+.1*std(delta_score{i}(:,j));
        [~,score_pks{i,j}]=findpeaks(delta_score{i}(:,j),'minpeakheight',threshold);
    end
end

%%

% express as a fraction of rando?

max_r=cell(1,length(use_photometry));
mean_r=cell(1,length(use_photometry));

for i=1:length(use_photometry)
    mat_photometry=markolab_vec2mat(zscore(use_photometry{i}.kin.ref.data),win_size,0);
    nwins=size(mat_photometry,2);
    max_r{i}=zeros(nwins,length(delta_win));
    for j=1:length(delta_win)
        mat_score=markolab_vec2mat(zscore(delta_score{i}(:,j)),win_size,0);
        for k=1:nwins
            r=xcorr(mat_photometry(:,k),mat_score(:,k),20,'coeff');
            max_r{i}(k,j)=max(r);
            mean_r{i}(k,j)=mean(r);
        end
    end
end

%%

upd=kinect_proctimer(nrands);
rnd_summary=[];

rnd_max_r=cell(1,length(use_photometry));
rnd_mean_r=cell(1,length(use_photometry));

for ii=1:nrands
    
    
    for i=1:length(use_photometry)
        
        mat_photometry=markolab_vec2mat(zscore(phase_rnds{i}.ref(:,ii)),win_size,0);
        nwins=size(mat_photometry,2);
        
        for j=1:length(delta_win)
            mat_score=markolab_vec2mat(zscore(delta_score{i}(:,j)),win_size,0);
            
            for k=1:nwins
                r=xcorr(mat_photometry(:,k),mat_score(:,k),20,'coeff');
                rnd_max_r{i}(k,j,ii)=max(r);
                rnd_mean_r{i}(k,j,ii)=mean(r);
            end
        end
    end
    
    upd(ii);
    
end

%%

% for plotting interpolate data, use gramm for the rest

save('timescale_analysis.mat','max_r','mean_r','score_pks','delta_score','delta_win','smooth_sig','delta_thresh','rnd_max_r','rnd_mean_r','-v7.3');

