% assumes we have data from master script

nboots=5e2;
maxlags=50;
delta_thresh=.05;
delta_wins=[1:2:40];
smooth_sig=2;
jl_bound_eps=.25;


%%

recon=U(:,1:rankcut)*S(1:rankcut,1:rankcut)*V(:,1:rankcut)'+mu;
recon=recon';
[gcamp_data,gcamp_ts]=kinect_analysis_proc_photometry(photometry_traces(1,:)+photometry_traces(2,:)*1e2,photometry_ts,'dff',1,'smooth_tau',.2);
gcamp_interp=interp1(gcamp_ts,gcamp_data,depth_ts,'linear','extrap');

%%

%n_components=round(4*log(size(recon,2))/(jl_bound_eps^2/2-jl_bound_eps^3/3));
n_components=200;
rps=kinect_gaussproj(zscore(recon),n_components);
rps=zscore(rps');
h=normpdf(-10:10,0,smooth_sig); % gauss smooth


%%

delta_score={};
changepoints={};
bootvals_gcamp={};
obs_r_gcamp={};
lags={};

for i=1:length(delta_wins)
  
    i
    rps_deltas=markolab_deltacoef(rps,delta_wins(i)); % delta coefficients, lag of 4 frames
    score=sum(abs(rps_deltas)>delta_thresh); % binarize deltas
    delta_score{i}=conv(score,h,'same'); % convolve
    thresh=mean(delta_score{i})+.5*std(delta_score{i});
    [~,changepoints{i}]=findpeaks(delta_score{i},'minpeakheight',thresh); 
    [bootvals_gcamp{i},~,obs_r_gcamp{i}]=kinect_analysis_xcorr_bootstr(zscore(gcamp_interp),zscore(delta_score{i}),...
        'nboots',nboots,'maxlags',maxlags);
    
end

%%


figs.timescale_scan=figure();

bootmu=cellfun(@mean,bootvals_gcamp);
bootstd=cellfun(@std,bootvals_gcamp);
corrvals=cellfun(@max,obs_r_gcamp);
med_ibi=cellfun(@(x) median(diff(x))*30/1e3,changepoints);
set(figs.timescale_scan,'defaultAxesColorOrder',[1 0 0;1 1 1]);

yyaxis left;
plot((corrvals-bootmu)./bootstd,'r.-','markersize',15);
ylabel('Max r (Z, bootstrap)');
set(gca,'ticklength',[0 0]);
yyaxis right;
plot(med_ibi,'g.-','markersize',15,'color',[1 1 1]);
set(gca,'ticklength',[0 0]);
% grid minor;
% set(gca,'MinorGridLineStyle','-');
box off;
set(figs.timescale_scan,'position',[300 300 250 200],'inverthardcopy','off','paperpositionmode','auto');
xlim([0 length(delta_wins)+1]);



