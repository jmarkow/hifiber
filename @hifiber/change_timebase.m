function change_timebase(OBJ,NEW_TIMESTAMPS)
%
%
%
%

valid_types={'dff','raw','baseline','baseline_rem','reref','reference'};

% interpolate all traces using the new timestamps

for i=1:length(OBJ)
	for j=1:length(OBJ(i).traces)
		trace_types=fieldnames(OBJ(i).traces(j));
		for k=1:length(trace_types)
			if any(strcmp(trace_types{k},valid_types)) & ~isempty(OBJ(i).traces(j).(trace_types{k}))

				nans=isnan(OBJ(i).timestamps)|isnan(OBJ(i).traces(j).(trace_types{k}));

				use_ts=OBJ(i).timestamps;
				use_trace=OBJ(i).traces(j).(trace_types{k});

				use_ts(nans)=[];
				use_trace(nans)=[];

				new_trace=interp1(use_ts,use_trace,NEW_TIMESTAMPS);
				OBJ(i).traces(j).(trace_types{k})=new_trace;

			end
		end
	end

	OBJ(i).timestamps=NEW_TIMESTAMPS;
	OBJ(i).metadata.fs=round(1./nanmean(diff(NEW_TIMESTAMPS)));

end
