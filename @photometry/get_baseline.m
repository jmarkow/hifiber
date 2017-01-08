function get_baseline(OBJ)
%
% detrends the data using a sliding window stat
% first use vec2mat to reformat data, then apply function to each window

upd=photometry.proc_timer(length(OBJ.traces));

for i=1:length(OBJ.traces)

	% window the data

	use_data=OBJ.traces(i).raw;

	nwin=round(OBJ.options.photometry.baseline_win.*OBJ.metadata.fs);
	if mod(nwin,2)==0
		nwin=nwin+1;
	end
	padded_data=[use_data(floor(nwin/2):-1:1);use_data;use_data(end-floor(nwin/2-1):1:end)];
	proc_mat=vec2mat(padded_data,nwin,nwin-1);

	% apply the rolling statistic

	baseline=OBJ.options.photometry.baseline_fcn(proc_mat)';
	OBJ.traces(i).baseline=baseline(:);
	upd(i);

end
