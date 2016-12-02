function kinect_analysis_extract_photometry(DATA,CH,REF,varargin)
%
%
%
%
%

if nargin<3 | isempty(REF)
	REF=[];
end

if nargin<2 | isempty(CH)
	CH=1;
end

if nargin<1 | isempty(DATA)
	if exist('../nidaq.dat','file')
		DATA='../nidaq.dat';
	else
		[filename,pathname]=uigetfile('*.dat');
		DATA=fullfile(pathname,filename);
	end
end

[opts,~,opts_names]=kinect_analysis_get_defaults('common','photometry');

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	if any(strcmp(varargin{i},opts_names))
		opts.(varargin{i})=varargin{i+1};
	end
end

% load in the NI-saved data

[pathname,~,~]=fileparts(DATA);

% read in the kinect timestamps

kinect_ts=kinect_read_csv(fullfile(pathname,'depth_ts.txt'));
kinect_ts=kinect_ts(:,2);
mat=kinect_nidaq2mat(DATA)';
fprintf('Found %i channels in %s\n',size(mat,2)-1,DATA);

if opts.dc_offset_split
	fprintf('Adding DC offsets...\n');
	fprintf('DC offset scale: %g\n',opts.dc_offset_scale);
	for i=1:length(CH)
		fprintf('CH %i DC CH %i\n',CH(i),CH(i)+opts.dc_offset_ch);
		mat(:,CH(i))=mat(:,CH(i))+mat(:,CH(i)+opts.dc_offset_ch)*opts.dc_offset_scale;
	end
end

photometry.raw.timestamps=mat(:,end);
photometry.raw.data=mat(:,CH);

if opts.invert_sig
	fprintf('Inverting signal...\n');
	photometry.raw.data=-photometry.raw.data;
end

photometry.raw.units='volts';
photometry.raw.labels=CH;

[photometry.proc.data, photometry.proc.data_baseline, photometry.proc.timestamps]=...
	kinect_analysis_proc_photometry(photometry.raw.data,photometry.raw.timestamps);
photometry.proc.labels=CH;
[nsamples,nchannels]=size(photometry.proc.data);
photometry.proc.units='volts over baseline';

% if we have a reference, use it!

if ~isempty(REF)

	signal_ch=setdiff(CH,REF);
	signal_n=length(signal_ch);

	photometry.ref.timestamps=photometry.proc.timestamps;
	photometry.ref.data=zeros(nsamples,signal_n);
	photometry.ref.labels=signal_ch;

	for i=1:signal_n
		photometry.ref.data(:,i)=fluolab_rereference(...
			photometry.proc.data(:,photometry.proc.labels==signal_ch(i)),...
			photometry.proc.data(:,photometry.proc.labels==REF));
		photometry.ref.data(:,i)=1e2*photometry.ref.data(:,i)...
			./photometry.proc.data_baseline(:,photometry.proc.labels==signal_ch(i));
	end

	if opts.clip
		photometry.ref.data(photometry.ref.data<0)=0;
	end

	photometry.ref.units='dF/F_0 (percent)';

end

% interpolate to the kinect timestamps, use reference however we like...

data_types=fieldnames(photometry);
data_types(strcmp(data_types,'raw'))=[];

for i=1:length(data_types)
	photometry.kin.(data_types{i}).data=interp1(...
		photometry.(data_types{i}).timestamps,...
		photometry.(data_types{i}).data,...
		kinect_ts,'linear','extrap');
    if isfield(photometry.(data_types{i}),'data_baseline')
        photometry.kin.(data_types{i}).data_baseline=interp1(...
            photometry.(data_types{i}).timestamps,...
            photometry.(data_types{i}).data_baseline,...
            kinect_ts,'linear','extrap');  
    end
	photometry.kin.(data_types{i}).labels=photometry.(data_types{i}).labels;
	photometry.kin.(data_types{i}).timestamps=kinect_ts;
	photometry.kin.(data_types{i}).units=photometry.(data_types{i}).units;
end

data_types=fieldnames(photometry);

for i=1:length(data_types)
	photometry.(data_types{i})=orderfields(photometry.(data_types{i}));
end

photometry.parameters=opts;


save('photometry.mat','photometry');

% load in the timestamps
