function PARAMS=get_demod_reference(SIG,SAMPLING_RATE,REF_FS,REF)
% way simpler than other implementations I've seen out there, just take complex-valued
% fft, freq/phase is equivalent to a least squares fit of a sinusoid, so there you go
%
% unless our signals are seriously corrupted/overlapping, we shouldn't need much
% guidance to do this #blessed

if nargin<4
	REF=[];
end

if nargin<3
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
	[~,loc]=max(mag_fft)

	% frequency guess

	freq_init=fvec(loc)
	phase_init=mod(phase_fft(loc),2*pi)-pi;

else

	freq_init=REF_FS; % filter and use zero cross something something for phase?
	phase_init=0;

end

ts=[0:sig_len-1]/SAMPLING_RATE;
ts=ts(~isnan(SIG));

% use all this crap to derive a least squares estimate of the reference
% amp is rms obvi, baseline shift is mean obvi, freq is peak of fft OBVI, and phase shift is whatevs

init_params=[nanstd(SIG) freq_init 0 nanmean(SIG)]

fit= @(b) b(1).*(sin(2*pi*ts*b(2)+b(3)))+b(4)-SIG(~isnan(SIG))';

% works pretty well son, good job

% blessed our lord non-linear least squares, may you guide your herd to the global minimum

PARAMS=lsqnonlin(fit,init_params,[init_params(1)*.5 freq_init-5 -pi 0],[init_params(2)*2 freq_init+5 pi 1]);
