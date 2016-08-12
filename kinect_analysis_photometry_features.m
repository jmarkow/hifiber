%%% script to perform basic analysis w/ scalars
% assumes scalars and photometry already loaded in, in addition 'recon' for the changepoints analysis

[~,changepoint_score,changepoints]=kinect_analysis_changepoints(recon);

nrands=1e3;
lags=200; % in frames
deriv_win=5;
lags_vec=[-lags:lags];
pcs=1:20;

%% get derivates for downstream stuff

use_photometry=zscore(photometry.kin.ref.data(:,1));
nsamples=length(use_photometry);
sig_fft=fft(use_photometry);
sig_amp=abs(sig_fft);

scr_theta=angle(fft(rand(nsamples,nrands)));
scr_data=real(ifft(repmat(sig_amp,[1 nrands]).*exp(1j.*scr_theta)));

obs_r.pcs=nan(2*lags+1,length(pcs));
boot_r.pcs=nan(2*lags+1,length(pcs),nrands);
obs_r.pcs_dt=nan(size(obs_r.pcs));
boot_r.pcs_dt=nan(size(obs_r.pcs_dt));

for i=1:length(pcs)

	% for each scalar, get the real xcorr and the scrambled xcorr

	obs_r.pcs(:,i)=xcorr(use_photometry,features.scores(pcs(i),:),lags,'coeff');

	for j=1:nrands
		boot_r.pcs(:,i,j)=xcorr(scr_data(:,j),features.scores(pcs(i),:),lags,'coeff');
	end

end

pcs_dt=markolab_deltacoef(features.scores(pcs,:),deriv_win);

for i=1:length(pcs)

	% for each scalar, get the real xcorr and the scrambled xcorr

	obs_r.pcs_dt(:,i)=xcorr(use_photometry,pcs_dt(i,:),lags,'coeff');
	for j=1:nrands
		boot_r.pcs_dt(:,i,j)=xcorr(scr_data(:,j),pcs_dt(i,:),lags,'coeff');
	end

end

obs_r.changepoint_score=xcorr(use_photometry,changepoint_score,lags,'coeff');
boot_r.changepoint_score=nan(2*lags+1,nrands);

for i=1:nrands
	boot_r.changepoint_score(:,i)=xcorr(scr_data(:,i),changepoint_score,lags,'coeff');
end

mu=mean(boot_r.changepoint_score,2);
sig=std(boot_r.changepoint_score,[],2);
obs_z.changepoint_score=(obs_r.changepoint_score-mu)./sig;

mu=mean(boot_r.pcs,3);
sig=std(boot_r.pcs,[],3);
obs_z.pcs=(obs_r.pcs-mu)./sig;

mu=mean(boot_r.pcs_dt,3);
sig=std(boot_r.pcs_dt,[],3);
obs_z.pcs_dt=(obs_r.pcs_dt-mu)./sig;

% format nice 'n purty

obs_r=orderfields(obs_r);
boot_r=orderfields(boot_r);
obs_z=orderfields(obs_z);

% do the same for smoothed deriv
