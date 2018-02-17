classdef hifiber < handle & matlab.mixin.SetGet
	% This class implements everything you need to extract raw data
	% collected from the Kinect v2 (with or without cable)

	% the essentials

	properties

		user_data

	end

	% stuff the user can see but can't modify without using a class method

	properties (GetAccess=public,SetAccess={?hifiber,?phanalysis,?kinect_extract})

		options
		timestamps
		traces
		references
		metadata
		status

	end

	% the completely hidden stuff

	properties (Access=private)

	end

	% have the pca object be a constant object, instance should be shared
	% across an array of objects...

	methods

		function obj=hifiber(DATA,TIMESTAMPS,FS,CHANNEL_NAMES)


			if nargin<4
				CHANNEL_NAMES={};
			end

			if nargin<3
				FS=[];
			end

			if nargin<2
				return;
			end

			% construct the object, then need methods for filtering, re-referencing
			% and matching other timestamps (e.g. from camera)

			% TODO:  status, steps, and automating everything similar to kinect extract

			if isvector(DATA)
				DATA=DATA(:);
			end

			obj.use_defaults;

			TIMESTAMPS=TIMESTAMPS(:);
			if isempty(FS)
				fs=1./nanmean(diff(TIMESTAMPS));
			else
				fs=FS;
			end

			[nsamples,nchannels]=size(DATA);

			obj.timestamps=TIMESTAMPS;
			obj.metadata.fs=round(fs);

			for i=1:nchannels
				obj.traces(i).raw=DATA(:,i);
				obj.traces(i).baseline=[];
				obj.traces(i).baseline_rem=[];
				obj.traces(i).dff=[];
        obj.traces(i).reref=[];
				obj.traces(i).name=sprintf('Ch%i',i);
				obj.traces(i).mod_freq=[];
				obj.traces(i).reference_channel=[];
			end

			obj.set_channel_names(CHANNEL_NAMES);
			obj.status.invert=false;

		end

		function s = saveobj(obj)
			use_names=properties(obj);
			for i=1:length(use_names)
				s.(use_names{i})=obj.(use_names{i});
			end
		end

	end

	methods(Static)

		% doesn't require the kinect_extract object
		upd=proc_timer(nloops,varargin)
		[new_ref,fun]=get_demod_reference(signal, ts, sampling_rate, ref_fs, ref_sig)
		[new_sync,fun]=get_sync_signal(signal, ts, sampling_rate)
		filt_data=bandpass(signal, fc, bw, sampling_rate, order, type, ripple, attenuation)
		mat=vec2mat(vec,nwin,noverlap)
		[proj,residuals,rescaled]=vector_rejection(sig,reference)

		function obj = loadobj(s)
				if isstruct(s)
					use_names=fieldnames(s);
					obj=hifiber;
					for i=1:length(use_names)
						obj.(use_names{i})=s.(use_names{i});
					end
				else
					obj=s;
				end

			end
	end

end
