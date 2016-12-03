%
%
%
%


win_size=900;
delta_win=10;
smooth_sig=.25;
delta_thresh=.1;
max_lag=90;
bin_smooth=2;
nrands=1e3;
ca_type_fields={'ref'};
dir_list={};

%%

%use_session=setdiff(8:22,[16 18:19]);
%use_session=setdiff(6:13,11);

% merge all of the usable data here???

load('arhmm_labels.mat','state_labels');
load('experiment_data_scores.mat','frame_idx','smooth_scores');
load('experiment_data_neural.mat','photometry');
load('experiment_data_pca.mat','features');
load('use_session.mat','use_session');

use_photometry=photometry(use_session);
use_labels=cell(1,length(state_labels));
use_scores=cell(1,length(use_photometry));
% now what!?!?!?!

%%

for i=1:length(use_photometry)
    use_labels{i}=state_labels{use_session(i)}(~isnan(frame_idx{use_session(i)}));
    use_scores{i}=smooth_scores{use_session(i)}(:,~isnan(frame_idx{use_session(i)}));
end

%%

% randomization

% phase randomization

if exist('phase_randomized_photometry.mat','file')~=2
    fprintf('Phase randomizing...');
    upd=kinect_proctimer(length(use_photometry)*length(ca_type_fields));
    counter=1;
    phase_rnds=cell(1,length(use_photometry));
    
    for i=1:length(use_photometry)
        
        ca_type_rnd.ref=use_photometry{i}.kin.ref.data(:,1);
        
        for j=1:length(ca_type_fields)
            
            tmp_fft=fft(ca_type_rnd.(ca_type_fields{j}));
            mag=abs(tmp_fft);
            ang=angle(tmp_fft);
            mag=mag(:);
            ang=ang(:);
            
            % generate nrands permutations
            
            [~,permidx]=sort(rand(numel(ang),nrands));
            
            % synthesize the random signals
            
            phase_rnds{i}.(ca_type_fields{j})=single(real(ifft(repmat(mag,[1 nrands]).*exp(1j.*ang(permidx)))));
            
            upd(counter);
            counter=counter+1;
            
        end
    end
    save('phase_randomized_photometry.mat','phase_rnds','-v7.3');
else
    fprintf('Loading phase randomizations...\n');
    load('phase_randomized_photometry.mat','phase_rnds');
end



%%

nsamples=numel(use_photometry{1}.kin.proc.data(:,1));
nwins=floor(nsamples/win_size);

% lag_mat.cpoint.gcamp=zeros(length(use_photometry)*nwins,max_lag*2+1);
% lag_mat.cpoint.autofluo=zeros(length(use_photometry)*nwins,max_lag*2+1);
% lag_mat.cpoint.ref=zeros(length(use_photometry)*nwins,max_lag*2+1);
counter=1;
npcs=size(smooth_scores{1},1);
lag_mat=[];

for i=1:npcs
    for j=1:length(ca_type_fields)
        lag_mat(i).(ca_type_fields{j})=zeros(length(use_photometry)*nwins,max_lag*2+1);
    end
end

for i=1:length(use_photometry)
    
    %     ca_type.gcamp=markolab_vec2mat(zscore(use_photometry{i}.kin.proc.data(:,1)./use_photometry{i}.kin.proc.data_baseline(:,1)),win_size,0);
    %     ca_type.autofluo=markolab_vec2mat(zscore(use_photometry{i}.kin.proc.data(:,2)./use_photometry{i}.kin.proc.data_baseline(:,2)),win_size,0);
    ca_type.ref=markolab_vec2mat(zscore(use_photometry{i}.kin.ref.data(:,1)),win_size,0);
    
    for j=1:npcs
        tmp_pc=markolab_vec2mat(zscore(use_scores{i}(j,:)),win_size,0);
        for k=1:length(ca_type_fields)
            for l=1:nwins
                lag_mat(j).(ca_type_fields{k})(nwins*(i-1)+l,:)=...
                    xcorr(ca_type.(ca_type_fields{k})(:,l),tmp_pc(:,l),max_lag,'coeff');
            end
        end
    end
end

%%

% phase rnd control

upd=kinect_proctimer(nrands);

for ii=1:nrands
    
    for i=1:length(beh_type_fields)
        for j=1:length(ca_type_fields)
            rnd_mat.(beh_type_fields{i}).(ca_type_fields{j})=zeros(length(use_photometry)*nwins,max_lag*2+1);
        end
    end
    
    for i=1:length(use_photometry)
        
        %         ca_type.gcamp=markolab_vec2mat(zscore(phase_rnds{i}.gcamp(:,ii)),win_size,0);
        %         ca_type.autofluo=markolab_vec2mat(zscore(phase_rnds{i}.autofluo(:,ii)),win_size,0);
        ca_type.ref=markolab_vec2mat(zscore(phase_rnds{i}.ref(:,ii)),win_size,0);
        
        beh_type.delta=markolab_vec2mat(zscore(delta_score{i}),win_size,0);
        beh_type.delta_bin=markolab_vec2mat(zscore(delta_score_bin{i}),win_size,0);
        beh_type.delta_bin_onsets=markolab_vec2mat(zscore(delta_score_bin_onsets{i}),win_size,0);
        beh_type.delta_bin_offsets=markolab_vec2mat(zscore(delta_score_bin_offsets{i}),win_size,0);
        beh_type.delta_bin_midpoints=markolab_vec2mat(zscore(delta_score_bin_midpoints{i}),win_size,0);
        beh_type.delta_model=markolab_vec2mat(zscore(model_changepoints{i}),win_size,0);
        
        for j=1:length(beh_type_fields)
            for k=1:length(ca_type_fields)
                for l=1:nwins
                    rnd_mat.(beh_type_fields{j}).(ca_type_fields{k})(nwins*(i-1)+l,:)=...
                        xcorr(ca_type.(ca_type_fields{k})(:,l),beh_type.(beh_type_fields{j})(:,l),max_lag,'coeff');
                end
            end
        end
    end
    
    % summarize randomization
    
    for i=1:length(beh_type_fields)
        for j=1:length(ca_type_fields)
            rnd_summary.(beh_type_fields{i}).(ca_type_fields{j})(ii,:)=mean(rnd_mat.(beh_type_fields{i}).(ca_type_fields{j}));
        end
    end
    
    upd(ii);
    
end

%%

save('changepoint_analysis.mat','lag_mat','rnd_summary');

