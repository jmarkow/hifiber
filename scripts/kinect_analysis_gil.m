%%% analysis of model output from Gil's experiment

%%

% load in the data

model_output=load('~/Desktop/test_model_output.mat');
scalars=load('~/Desktop/test_data.mat');
scalars=orderfields(scalars);
states=unique(model_output.labels);
states(isnan(states))=[];
scalars.data=scalars.data';
max_plot=100;
rotate_plot=true;
surface_plot=false;
show_endpoints=false;
smoothing=8;
duration_min=5;
colors=[1 0 0;0 0 0;0 1 0];
colors_hsv=rgb2hsv(colors);

%%

mouse_names=cellstr(model_output.mouse_names);
dorsal_stim=~cellfun(@isempty,strfind(mouse_names,'Dorsal'));
ventral_stim=~cellfun(@isempty,strfind(mouse_names,'Ventral'));
no_stim=~cellfun(@isempty,strfind(mouse_names,'noStim'));
states=[4 5];

%%

for ii=1:length(states)

	syllable_id=states(ii); % so far 7 and 12 look good

	hits=model_output.labels==syllable_id;

	if mean(hits)<.01
		continue;
	end

	idx=1:length(model_output.labels)-1;
	rising_edges=find(~hits(idx)&hits(idx+1));
	falling_edges=find(hits(idx)&~hits(idx+1));

	if length(rising_edges)~=length(falling_edges)
		continue;
	end

	use_idx=[rising_edges(:) falling_edges(:)];
	durations=diff(use_idx,[],2);

	include_fields={'angle','head_x_mm','head_y_mm','tail_x_mm','tail_y_mm','centroid_x_mm','centroid_y_mm','data'};
	examples=struct();

	for i=1:size(use_idx,1)
		if ~isempty(strfind(mouse_names{use_idx(i,1)},'Dorsal'))
			condition='dorsal';
		elseif ~isempty(strfind(mouse_names{use_idx(i,1)},'Ventral'))
			condition='ventral';
		else
			condition='nostim';
		end
		for j=1:length(include_fields)
			examples(i).(include_fields{j})=scalars.(include_fields{j})(:,use_idx(i,1):use_idx(i,2));
		end
		examples(i).condition=condition;
		examples(i).duration=size(examples(i).data,2);
	end

	use_examples=examples(durations>=duration_min);
	use_conditions={use_examples(:).condition};

	% get the n of each condition

	idx_dorsal=find(strcmp(use_conditions,'dorsal'));
	idx_nostim=find(strcmp(use_conditions,'nostim'));
	idx_ventral=find(strcmp(use_conditions,'ventral'));

	if length(idx_dorsal)>max_plot
		idx_dorsal=randsample(idx_dorsal,max_plot,'false');
	end

	if length(idx_nostim)>max_plot
		idx_nostim=randsample(idx_nostim,max_plot,'false');
	end

	if length(idx_ventral)>max_plot
		idx_ventral=randsample(idx_ventral,max_plot,'false');
	end

	if length(idx_dorsal)<1 | length(idx_nostim)<1 | length(idx_ventral)<1
		continue;
	end

	use_examples=use_examples([idx_dorsal(:);idx_nostim(:);idx_ventral(:)]);
	fig=[];

	fig(ii)=figure('name',sprintf('Syllable %i',states(ii)));
	plot_variables={'tail_x_mm','tail_y_mm','centroid_x_mm','centroid_y_mm','head_x_mm','head_y_mm'};

	dorsal_counter=1;
	ventral_counter=1;
	nostim_counter=1;

	dorsal_colors=interp1([1 length(idx_dorsal)/2 length(idx_dorsal)],colors_hsv,1:length(idx_dorsal));
	dorsal_colors=hsv2rgb(dorsal_colors);

	ventral_colors=interp1([1 length(idx_ventral)/2 length(idx_ventral)],colors_hsv,1:length(idx_ventral));
	ventral_colors=hsv2rgb(ventral_colors);

	nostim_colors=interp1([1 length(idx_nostim)/2 length(idx_nostim)],colors_hsv,1:length(idx_nostim));
	nostim_colors=hsv2rgb(nostim_colors);

	max_duration=prctile((cat(1,use_examples(:).duration)),95);

	%dorsal_colors=colormap([colors '(' num2str(length(idx_dorsal)) ')']);
	%ventral_colors=colormap([colors '(' num2str(length(idx_ventral)) ')']);
	%nostim_colors=colormap([colors '(' num2str(length(idx_nostim)) ')']);

	for i=1:length(use_examples)

		if strcmp(use_examples(i).condition,'dorsal')
			linecolor=dorsal_colors(dorsal_counter,:);
			dorsal_counter=dorsal_counter+1;
		elseif strcmp(use_examples(i).condition,'ventral')
			linecolor=ventral_colors(ventral_counter,:);
			ventral_counter=ventral_counter+1;
		else
			linecolor=nostim_colors(nostim_counter,:);
			nostim_counter=nostim_counter+1;
		end

		for j=1:2:length(plot_variables)

			nsamples=size(use_examples(i).data,2);
			plot_centroid=[use_examples(i).(plot_variables{j});use_examples(i).(plot_variables{j+1})];
			plot_centroid=plot_centroid-repmat(plot_centroid(:,1),[1 nsamples]);

			if smoothing>0
				plot_centroid=filter(ones(smoothing,1)/smoothing,1,plot_centroid')';
			end

			if rotate_plot
				% rotate so the animal starts facing to the right
				angrad=-(mod(use_examples(i).angle,2*pi));
				%angrad=use_examples(i).angle;
				rotation=[cos(angrad(1)) -sin(angrad(1)); sin(angrad(1)) cos(angrad(1))];
				for k=1:nsamples
					plot_centroid(:,k)=rotation*plot_centroid(:,k);
				end
			end

			% don't need to flip x axis (positive changes right, negative left already)

			x=plot_centroid(1,:);

			% flip y axis (positive changes are down in movies, negative up)

			y=-plot_centroid(2,:);
			z=zeros(size(x));
			col=min([1:size(x,2)]./max_duration,1);

			if strcmp(use_examples(i).condition,'dorsal')
				subplot(3,3,floor(j/2)+1);
			elseif strcmp(use_examples(i).condition,'ventral')
				subplot(3,3,3+floor(j/2)+1);
			else
				subplot(3,3,6+floor(j/2)+1);
			end

			if surface_plot
				surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',.5);
				colormap(jet);
			else
				plot(x,y,'k-','linewidth',.5,'color',linecolor);
			end

			hold on;

			if show_endpoints
				scatter(x(1),y(1),5,'bo','filled','markeredgecolor','k',...
					'markerfacecolor',linecolor);
				scatter(x(end),y(end),5,'ro','filled','markeredgecolor','k',...
					'markerfacecolor',linecolor);
			end
		end

	end

	ax=[];

	xlimits=[inf -inf];
	ylimits=[inf -inf];

	for i=1:9

		ax(i)=subplot(3,3,i);
		axis square
		axis off

		xlims=xlim();
		ylims=ylim();

		if xlims(1)<xlimits(1), xlimits(1)=xlims(1); end
		if xlims(2)>xlimits(2), xlimits(2)=xlims(2); end
		if ylims(1)<ylimits(1), ylimits(1)=ylims(1); end
		if ylims(2)>ylimits(2), ylimits(2)=ylims(2); end

	end

	xlimit=round(max(abs(xlimits))/10)*10;
	ylimit=round(max(abs(ylimits))/10)*10;

	xlimit=max(xlimit,ylimit);
	ylimit=max(xlimit,ylimit);

	xlimit=150;
	ylimit=150;

	linkaxes(ax,'xy');
	axis([-xlimit xlimit -ylimit ylimit]);
	ticklength=7;
	offset=2;
	spacing=50;

	for i=1:9

		subplot(3,3,i);
		h=[];
		h(1)=line([-spacing spacing],[-ylimit-offset -ylimit-offset]);
		h(2)=line([-(xlimit+offset) -(xlimit+offset)],[-spacing spacing]);
		h(3)=line([-spacing -spacing],[-ylimit-offset -ylimit-(offset+ticklength)]);
		h(4)=line([0 0],[-ylimit-offset -ylimit-(offset+ticklength)]);
		h(5)=line([spacing spacing],[-ylimit-offset -ylimit-(offset+ticklength)]);
		h(6)=line([-(xlimit+(offset+ticklength)) -(xlimit+offset)],[-spacing -spacing]);
		h(7)=line([-(xlimit+(offset+ticklength)) -(xlimit+offset)],[0 0]);
		h(8)=line([-(xlimit+(offset+ticklength)) -(xlimit+offset)],[spacing spacing]);

		for j=1:length(h)
			set(h(j),'clipping','off','color','k');
		end

	end

	%axis([-200 200 -200 200]);
	% subplot(3,3,9);
	% h=line([xlimit-25 xlimit],[-ylimit -ylimit],'color','k');
	% h2=line([xlimit xlimit],[-ylimit -ylimit+25],'color','k');
	% set(h,'clipping','off');
	% set(h2,'clipping','off');

	subplot(3,3,1);
	t=text(-xlimit-xlimit/2,-ylimit/2,[ 'Dorsal (' num2str(dorsal_counter-1) ')']);
	set(t,'Rotation',90);
	t=text(-xlimit/4,ylimit+ylimit/10,'Tail');

	subplot(3,3,2);
	t=text(-xlimit/4,ylimit+ylimit/10,'Centroid');
	subplot(3,3,3);
	t=text(-xlimit/4,ylimit+ylimit/10,'Head');

	subplot(3,3,4);
	t=text(-xlimit-xlimit/2,-ylimit/2,[ 'Ventral (' num2str(ventral_counter-1) ')' ]);
	set(t,'Rotation',90);

	subplot(3,3,7);
	t=text(-xlimit-xlimit/2,-ylimit/2,[ 'No Stim (' num2str(nostim_counter-1) ')' ]);
	set(t,'Rotation',90);
	markolab_multi_fig_save(fig(ii),'~/Desktop/for_gil/',...
		sprintf('syllable_%i_rotation_%i_endpoints_%i_nsamples_%i_surface_%i',...
			states(ii),rotate_plot,show_endpoints,max_plot,surface_plot),'eps,png,fig','renderer','painters');

end
