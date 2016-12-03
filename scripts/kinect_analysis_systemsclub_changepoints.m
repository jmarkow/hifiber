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
beh_type_fields={'delta','delta_bin','delta_model','delta_bin_onsets','delta_bin_offsets','delta_bin_midpoints'};
ca_type_fields={'ref','proc1','proc2'};
dir_list={};

%%

%use_session=setdiff(8:22,[16 18:19]);
%use_session=setdiff(6:13,11);

% merge all of the usable data here???

load('arhmm_labels.mat','state_labels');
load('experiment_data_scores.mat','frame_idx');
load('experiment_data_neural.mat','photometry');
load('experiment_data_rps.mat','rps');
load('use_session.mat','use_session');
use_photometry=photometry(use_session);
use_rps=rps(use_session);
use_labels=state_labels(use_session);
use_frame_idx=frame_idx(use_session);

if isfield(use_photometry{i}.kin,'ref')
    ca_type_fields={'ref'};
else
    ca_type_fields={'proc1','proc2'};
end

% now what!?!?!?!

%% first get changepoints, doo dah

if exist('delta_score.mat','file')~=2
    timer_upd=kinect_proctimer(length(delta_win)*length(use_photometry));
    delta_score=cell(1,length(use_photometry));
    
    counter=1;
    for i=1:length(use_photometry)
        
        delta_score{i}=zeros(size(use_rps{i},1),length(delta_win));
        tmp_rps=zscore(zscore(use_rps{i})');
        %tmp_rps=hampel(tmp_rps',5,3)';
        deltas=markolab_deltacoef(tmp_rps,delta_win); % delta coefficients, lag of 4 frames
        delta_score{i}=sum(abs(deltas)>delta_thresh); % binarize deltas
        h=normpdf([-10:10],0,smooth_sig); % gauss smooth
        delta_score{i}=conv(delta_score{i},h,'same'); % convolve
        timer_upd(counter);
        counter=counter+1;
        
    end
    save('delta_score.mat','delta_score');
else
    fprintf('Loading delta scores...\n');
    load('delta_score.mat','delta_score');
end

%%

delta_score_bin=cell(size(delta_score));
delta_score_bin_onsets=cell(size(delta_score));
delta_score_bin_offsets=cell(size(delta_score));
delta_score_bin_midpoints=cell(size(delta_score));

smooth_kernel=normpdf([-10:10],0,bin_smooth);

for i=1:length(use_photometry)
    tmp=delta_score{i};
    thresh=mean(tmp)+.1*std(tmp);
    delta_score_bin{i}=zeros(size(tmp));
    delta_score_bin_onsets{i}=zeros(size(tmp));
    delta_score_bin_offsets{i}=zeros(size(tmp));
    delta_score_bin_midpoints{i}=zeros(size(tmp));
    
    [~,locs]=findpeaks(tmp,'minpeakheight',thresh);
    %
    %     % get rise times for boundaries
    %     idx=1:length(tmp)-1;
    %     rise_times=tmp(idx)<thresh&tmp(idx+1)>thresh;
    %
    %delta_score_bin{i}=rise_times;
    locs=unique([1;locs(:);length(tmp)]);
    onsets=locs(1:2:length(locs));
    offsets=locs(2:2:length(locs));
    
    if onsets(end)==length(tmp)
        onsets(end)=onsets(end)-1;
        offsets(end+1)=length(tmp);
    end
    
    delta_score_bin{i}(locs)=1;
    delta_score_bin{i}=conv(delta_score_bin{i},smooth_kernel,'same');
    delta_score_bin_onsets{i}(onsets)=1;
    delta_score_bin_onsets{i}=conv(delta_score_bin_onsets{i},smooth_kernel,'same');
    delta_score_bin_offsets{i}(offsets)=1;
    delta_score_bin_offsets{i}=conv(delta_score_bin_offsets{i},smooth_kernel,'same');
    delta_score_bin_midpoints{i}(round(median([onsets(:)';offsets(:)'])))=1;
    delta_score_bin_midpoints{i}=conv(delta_score_bin_midpoints{i},smooth_kernel,'same');
end


model_changepoints=cell(size(delta_score));

% first map states to original vectors then

for i=1:length(delta_score)
    use_labels{i}=state_labels{use_session(i)}(~isnan(frame_idx{use_session(i)}));
    model_changepoints{i}=single([0 abs(diff(use_labels{i}))>0]);
    model_changepoints{i}=conv(model_changepoints{i},smooth_kernel,'same');
end


%%

% randomization

% phase randomization

if exist('phase_randomized_photometry.mat','file')~=2
    fprintf('Phase randomizing...\n');
    upd=kinect_proctimer(length(use_photometry)*length(ca_type_fields));
    counter=1;
    phase_rnds=cell(1,length(use_photometry));
    
    for i=1:length(use_photometry)
        
        
        if isfield(use_photometry{i}.kin,'ref')
            ca_type_rnd.ref=use_photometry{i}.kin.ref.data(:,1);
        else
            ca_type_rnd.proc1=use_photometry{i}.kin.proc.data(:,1);
            ca_type_rnd.proc2=use_photometry{i}.kin.proc.data(:,2);
        end
        
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


for i=1:length(beh_type_fields)
    for j=1:length(ca_type_fields)
        lag_mat.(beh_type_fields{i}).(ca_type_fields{j})=zeros(length(use_photometry)*nwins,max_lag*2+1);
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
    
    beh_type.delta=markolab_vec2mat(zscore(delta_score{i}),win_size,0);
    beh_type.delta_bin=markolab_vec2mat(zscore(delta_score_bin{i}),win_size,0);
    beh_type.delta_bin_onsets=markolab_vec2mat(zscore(delta_score_bin_onsets{i}),win_size,0);
    beh_type.delta_bin_offsets=markolab_vec2mat(zscore(delta_score_bin_offsets{i}),win_size,0);
    beh_type.delta_bin_midpoints=markolab_vec2mat(zscore(delta_score_bin_midpoints{i}),win_size,0);
    beh_type.delta_model=markolab_vec2mat(zscore(model_changepoints{i}),win_size,0);
    
    for j=1:length(beh_type_fields)
        for k=1:length(ca_type_fields)
            for l=1:nwins
                lag_mat.(beh_type_fields{j}).(ca_type_fields{k})(nwins*(i-1)+l,:)=...
                    xcorr(ca_type.(ca_type_fields{k})(:,l),beh_type.(beh_type_fields{j})(:,l),max_lag,'coeff');
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
    
    for i=1:length(beh_type_fields)
        for j=1:length(ca_type_fields)
            rnd_mat.(beh_type_fields{i}).(ca_type_fields{j})=zeros(length(use_photometry)*nwins,max_lag*2+1);
        end
    end
    
    for i=1:length(use_photometry)
        
        %         ca_type.gcamp=markolab_vec2mat(zscore(phase_rnds{i}.gcamp(:,ii)),win_size,0);
        %         ca_type.autofluo=markolab_vec2mat(zscore(phase_rnds{i}.autofluo(:,ii)),win_size,0);
        
        if isfield(use_photometry{i}.kin,'ref')
            ca_type.ref=markolab_vec2mat(zscore(phase_rnds{i}.ref(:,ii)),win_size,0);
        else
            ca_type.proc1=markolab_vec2mat(zscore(phase_rnds{i}.proc1(:,ii)),win_size,0);
            ca_type.proc2=markolab_vec2mat(zscore(phase_rnds{i}.proc2(:,ii)),win_size,0);
        end
%         ca_type.ref=markolab_vec2mat(zscore(phase_rnds{i}.ref(:,ii)),win_size,0);

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

