function demodulate(OBJ)
%
%
%

% use supplied references to demodulate each channel using lock-in principles

% (i.e. multiply by ref and 90 deg phase-shift)

% in other words, hypot it out and smooth yeahhhhhh
fprintf('Demodulating signals...\n');

for i=1:length(OBJ)

	demod_samples=round(OBJ(i).metadata.fs*OBJ(i).options.demod_tau);
	demod_kernel=ones(demod_samples,1)/demod_samples;

	ntraces=length(OBJ(i).traces);
	traces=cat(2,OBJ(i).traces(:).raw);
	mod_freq=cat(2,OBJ(i).traces(:).mod_freq);

	OBJ(i).traces=[];
	OBJ(i).traces=[];

	counter=1;
	upd=hifiber.proc_timer(ntraces*length(OBJ(i).references));

	% just set traces to empty while we do this perhaps...
	for j=1:ntraces

		% remove traces, and re-populate with the demodded version
		% do all to all, and bandpass as needed...

		for k=1:length(OBJ(i).references)

			% bandpass for the reference fs

			use_data=hifiber.bandpass(traces(:,j),mod_freq(k),...
				OBJ(i).options.mod_bandpass_bw,OBJ(i).metadata.fs);

			% TODO:  much more intelligent handling of the pads fool

			mult_x=use_data(:).*OBJ(i).references(k).x(:);
			mult_y=use_data(:).*OBJ(i).references(k).y(:);

			% linear conv zero pads the end of u, be sure to nan that ish out

			prod_x=conv(mult_x,demod_kernel,'same');
			prod_y=conv(mult_y,demod_kernel,'same');

			OBJ(i).traces(counter).raw=hypot(prod_x,prod_y);
			OBJ(i).traces(counter).raw(1:(demod_samples+1))=nan;
			OBJ(i).traces(counter).raw(end-(demod_samples+1):end)=nan;

			OBJ(i).traces(counter).type='Data';
			OBJ(i).traces(counter).name=sprintf('Ch %i demod ref %i',j,k);

			upd(counter);
			counter=counter+1;

		end


	end
end
