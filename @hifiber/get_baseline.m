function get_baseline(OBJ)
%
% detrends the data using a sliding window stat
% first use vec2mat to reformat data, then apply function to each window


for i=1:length(OBJ)

	if length(OBJ(i).traces)>0
		upd=hifiber.proc_timer(length(OBJ(i).traces));
	end

	for j=1:length(OBJ(i).traces)

		% window the data

		use_data=OBJ(i).traces(j).raw;

		nwin=round(OBJ(i).options.baseline_win.*OBJ(i).metadata.fs);
		if mod(nwin,2)==0
			nwin=nwin+1;
		end

		%padded_data=[use_data(floor(nwin/2):-1:1);use_data;use_data(end-floor(nwin/2-1):1:end)];
		padded_data=[nan(floor(nwin/2),1)*use_data(1);use_data;nan(floor(nwin/2),1)*use_data(end)];
		proc_mat=hifiber.vec2mat(padded_data,nwin,nwin-1);

		% apply the rolling statistic

		baseline=OBJ(i).options.baseline_fcn(proc_mat)';

		if OBJ(i).options.baseline_post_smooth>0
			baseline=smooth(baseline,OBJ(i).options.baseline_post_smooth*OBJ(i).metadata.fs);
		end

		OBJ(i).traces(j).baseline=baseline(:);
		upd(j);

	end
end
