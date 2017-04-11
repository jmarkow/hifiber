function change_timebase(OBJ,NEW_TIMESTAMPS)
%
%
%
%

% interpolate all traces using the new timestamps

for i=1:length(OBJ.traces)
	trace_types=fieldnames(OBJ.traces(i));
	for j=1:length(trace_types)
		new_trace=interp1(OBJ.timestamps,OBJ.traces(i).(trace_types{j}),NEW_TIMESTAMPS,'spline');
		OBJ.traces(i).(trace_types{j})=new_trace;
	end
end

OBJ.timestamps=NEW_TIMESTAMPS;
OBJ.metadata.fs=round(1./mean(diff(NEW_TIMESTAMPS)));
