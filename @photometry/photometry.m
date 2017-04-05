classdef photometry < handle & matlab.mixin.SetGet
	% This class implements everything you need to extract raw data
	% collected from the Kinect v2 (with or without cable)

	% the essentials

	properties

	end

	% stuff the user can see but can't modify without using a class method

	properties (GetAccess=public,SetAccess=private)

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

		function obj=photometry(DATA,TIMESTAMPS,CHANNEL_NAMES)

			if nargin<3
				CHANNEL_NAMES={};
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
			fs=1./mean(diff(TIMESTAMPS));

			[nsamples,nchannels]=size(DATA);

			obj.timestamps=TIMESTAMPS;
			obj.metadata.fs=round(fs);

			for i=1:nchannels
				obj.traces(i).raw=DATA(:,i);
				obj.metadata.traces(i).type='data';
				obj.metadata.traces(i).name=sprintf('Ch%i',i);
				obj.metadata.traces(i).mod_freq=[];
			end

			obj.set_channel_names(CHANNEL_NAMES);
			obj.status.invert=false;

		end

		% function saveobj(obj)
		% 	s.options=obj.options;
		% 	s.traces=obj.traces;
		% 	s.metadata=obj.metadata;
		% 	s.timestamps=obj.timestamps;
		% end

	end

	methods(Static)

		% doesn't require the kinect_extract object
		upd=proc_timer(nloops,varargin)
		[new_ref,fun]=get_demod_reference(signal, ts, sampling_rate, ref_fs, ref_sig)
		[new_sync,fun]=get_sync_signal(signal, ts, sampling_rate)
		filt_data=bandpass(signal, fc, bw, sampling_rate)
		% function obj=loadobj(s)
		% 	if isstruct(s)
		% 		newobj=photometry;
		% 		newobj.options=s.options;
		% 		newobj.traces=s.traces;
		% 		newobj.metadata=s.metadata;
		% 		newobj.timestamps=s.timestamps;
		% 		obj=newobj;
		% 	end
		%
		% end

	end

end
