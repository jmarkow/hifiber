%
%
%
%


win_size=900;
smooth_sig=.25;
delta_thresh=.1;
max_lag=90;
bin_smooth=2;
nrands=1e3;
dir_list={};

%%

%use_session=setdiff(8:22,[16 18:19]);
%use_session=setdiff(6:13,11);

% merge all of the usable data here???

load('experiment_data_scores.mat','frame_idx','smooth_scores');
load('experiment_data_neural.mat','photometry');
load('experiment_data_pca.mat','features');
load('use_session.mat','use_session');
load('phase_randomized_photometry.mat','phase_rnds');

use_photometry=photometry(use_session);
use_scores=cell(1,length(use_photometry));
% now what!?!?!?!

if isfield(use_photometry{1}.kin,'ref')
    ca_type_fields={'ref'};
else
    ca_type_fields={'proc1','proc2'};
end

%%

for i=1:length(use_photometry)
    use_scores{i}=smooth_scores{use_session(i)}(:,~isnan(frame_idx{use_session(i)}));
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
    
    if isfield(use_photometry{i}.kin,'ref')
        ca_type.ref=markolab_vec2mat(zscore(use_photometry{i}.kin.ref.data(:,1)),win_size,0);
    else
        ca_type.proc1=markolab_vec2mat(use_photometry{i}.kin.proc.data(:,1),win_size,0);
        ca_type.proc2=markolab_vec2mat(use_photometry{i}.kin.proc.data(:,2),win_size,0);
    end
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
rnd_summary=[];

for ii=1:nrands
    
    rnd_mat=[];
    
    for i=1:npcs
        for j=1:length(ca_type_fields)
            rnd_mat(i).(ca_type_fields{j})=zeros(length(use_photometry)*nwins,max_lag*2+1);
        end
    end
    
    for i=1:length(use_photometry)
        
        %         ca_type.gcamp=markolab_vec2mat(zscore(phase_rnds{i}.gcamp(:,ii)),win_size,0);
        %       ca_type.autofluo=markolab_vec2mat(zscore(phase_rnds{i}.autofluo(:,ii)),win_size,0);
        if isfield(use_photometry{i}.kin,'ref')
            ca_type.ref=markolab_vec2mat(zscore(phase_rnds{i}.ref(:,ii)),win_size,0);
        else
            ca_type.proc1=markolab_vec2mat(phase_rnds{i}.proc1(:,ii),win_size,0);
            ca_type.proc2=markolab_vec2mat(phase_rnds{i}.proc2(:,ii),win_size,0);
        end
         
        
        for j=1:npcs
            tmp_pc=markolab_vec2mat(zscore(use_scores{i}(j,:)),win_size,0);
            for k=1:length(ca_type_fields)
                for l=1:nwins
                    rnd_mat(j).(ca_type_fields{k})(nwins*(i-1)+l,:)=...
                        xcorr(ca_type.(ca_type_fields{k})(:,l),tmp_pc(:,l),max_lag,'coeff');
                end
            end
        end
    end
    
    % summarize randomization
    
    for i=1:npcs
        for j=1:length(ca_type_fields)
            rnd_summary(i).(ca_type_fields{j})(ii,:)=mean(rnd_mat(i).(ca_type_fields{j}));
        end
    end
    
    upd(ii);
    
end

%%

save('pca_analysis.mat','lag_mat','rnd_summary');

