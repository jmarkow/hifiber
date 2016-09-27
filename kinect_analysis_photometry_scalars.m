function [obs_r,obs_p,obs_z,boot_r]=kinect_analysis_photometry_scalars(DIR,varargin)
%%% script to perform basic analysis w/ scalars
% assumes scalars and photometry already loaded in

if nargin<1 | isempty(DIR)
	DIR=fullfile(pwd,'analysis');
end

load(fullfile(DIR,'features.mat'),'scalars');
load(fullfile(DIR,'photometry.mat'),'photometry');

% option to set random seed for reproducibility?

file_save=true;
nrands=1e4;
lags=100; % in frames
deriv_win=5;
lags_vec=[-lags:lags];

%% get derivates for downstream stuff

use_scalars=scalars;
scalar_names=fieldnames(use_scalars);
exclude=~cellfun(@isempty,regexp(scalar_names,'(centroid|skewness|theta)'));
scalar_names(exclude)=[];

for i=1:length(scalar_names)
	new_name=[ scalar_names{i} '_dt' ];
	scalar_names{end+1}=new_name;
	use_scalars.(new_name)=markolab_deltacoef(scalars.(scalar_names{i}),deriv_win)';
end

scalar_names(strcmp(scalar_names,'angle'))=[]; % will need to come back to do linear-circular corr

% phase scramble photometry signal for stats

fprintf('Phase scrambling photometry signal\n');

use_photometry=zscore(photometry.kin.ref.data(:,1));
nsamples=length(use_photometry);
sig_fft=fft(use_photometry);
sig_amp=abs(sig_fft);

scr_theta=angle(fft(rand(nsamples,nrands)));
scr_data=real(ifft(repmat(sig_amp,[1 nrands]).*exp(1j.*scr_theta)));

for i=1:length(scalar_names)

	% for each scalar, get the real xcorr and the scrambled xcorr

	fprintf('Analyzing feature %i of %i: %s\n',i,length(scalar_names),scalar_names{i});

	obs_r.(scalar_names{i})=xcorr(use_photometry,use_scalars.(scalar_names{i}),lags,'coeff');
	boot_r.(scalar_names{i})=nan(lags*2+1,nrands);

	for j=1:nrands
		boot_r.(scalar_names{i})(:,j)=...
			xcorr(scr_data(:,j),use_scalars.(scalar_names{i}),lags,'coeff');
	end

	mu=mean(boot_r.(scalar_names{i}),2);
	sig=std(boot_r.(scalar_names{i}),[],2);
	obs_z.(scalar_names{i})=(obs_r.(scalar_names{i})-mu)./sig;

	% right tail, alternative is that correlation is > bootstrap

	obs_p.right.(scalar_names{i})=mean(repmat(obs_r.(scalar_names{i}),[1 nrands])...
		<=repmat(max(boot_r.(scalar_names{i})),[2*lags+1 1]),2)+1/nrands;

	% left tail, alternative is that correlation is < bootstrap

	obs_p.left.(scalar_names{i})=mean(repmat(obs_r.(scalar_names{i}),[1 nrands])...
		>=repmat(min(boot_r.(scalar_names{i})),[2*lags+1 1]),2)+1/nrands;

end

% format nice 'n purty-like

obs_p.right=orderfields(obs_p.right);
obs_p.left=orderfields(obs_p.left);
obs_r=orderfields(obs_r);
boot_r=orderfields(boot_r);
obs_z=orderfields(obs_z);

if file_save
	save(fullfile(DIR,'analysis_scalars.mat'),'obs_p','obs_r','obs_z','boot_r');
end
