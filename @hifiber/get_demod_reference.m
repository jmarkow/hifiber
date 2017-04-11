function [PARAMS,FIT_FUN]=get_demod_reference(SIG,TS,SAMPLING_RATE,REF_FS,REF)
% way simpler than other implementations I've seen out there, just take complex-valued
% fft, freq/phase is equivalent to a least squares fit of a sinusoid, so there you go
%
% unless our signals are seriously corrupted/overlapping, we shouldn't need much
% guidance to do this #blessed

if nargin<5
	REF=[];
end

if nargin<5
	REF_FS=[];
end

% chop out nans and shit

SIG(isnan(SIG))=0;

% take an fft or whatever it's called

sig_len=length(SIG);

if isempty(REF_FS)

	use_sig=SIG;
	use_sig(isnan(use_sig))=0;

	comp_fft=fft(use_sig-mean(use_sig));

	% chop off the second half, rescale energy cuz you know you want to!

	comp_fft=comp_fft(1:sig_len/2+1);

	% conserve energy says perseval

	comp_fft(2:end-1)=2*comp_fft(2:end-1);

	% mag and phase bwooooomp

	mag_fft=abs(comp_fft/sig_len);
	mag_fft(1)=0;
	% fix ze angles d00d

	phase_fft=mod(unwrap(angle(comp_fft)),2*pi);

	% ref is now sin(2*pi*f0*t+phase_shift) and quadrature is sin(2*pi*f0*t+phase_shift+pi/2)

	fvec=SAMPLING_RATE*[0:(sig_len/2)]/sig_len;
	[~,loc]=max(mag_fft);

	% frequency guess

	freq_init=fvec(loc);
	%phase_init=phase_fft(loc);
	phase_init=0;

else

	freq_init=REF_FS; % filter and use zero cross something something for phase?
	phase_init=0;

end

% use all this crap to derive a least squares estimate of the reference
% amp is rms obvi, baseline shift is mean obvi, freq is peak of fft OBVI, and phase shift is whatevs

init_params=[nanstd(SIG) freq_init 0 nanmean(SIG)];

% exclude any bad datersss
TS=TS(~isnan(SIG));
SIG=SIG(~isnan(SIG));

TS=TS(:);
SIG=SIG(:);

fprintf('Initial parameters: %g\n',init_params)

FIT_FUN= @(b,x) b(1).*(sin(2*pi*x*b(2)+b(3)))+b(4);
obj_fun= @(b) FIT_FUN(b,TS)-SIG;

% works pretty well son, good job

% blessed our lord non-linear least squares, may you guide your herd to the global minimum
% get an initial phase shift by fitting the sinusoid and checking angle difference

init_fit_angle=angle(hilbert(FIT_FUN(init_params,TS)));
data_angle=angle(hilbert(SIG));
phase_diff=angle(mean(exp(1j.*(init_fit_angle(:)-data_angle(:)))));
init_params(3)=angle(exp(1j.*(-phase_diff))); % it's 0 to start, so just get angle diff from 0 dawg
opts=optimset('Display','off');
%fprintf('Detected a phase difference of %g(rads), correcting before optimization\n',phase_diff);
%opts=optimset('FinDiffRelStep',[.01 .05 .1 .2],'TolFun',1e-8,'TolX',1e-8);

PARAMS=lsqnonlin(obj_fun,init_params,[init_params(1)*.75 freq_init-20 -2*pi 0],[init_params(2)*1.25 freq_init+20 2*pi 1],opts);
