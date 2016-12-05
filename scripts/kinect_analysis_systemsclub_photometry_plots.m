%
%
%
%
%
%
%

% changepoints stuff

% get the std from the bootstrap

% maybe just insert a loop for this


camera_fs=30;
conditions=fieldnames(lag_mat);
conditions(strcmp(conditions,'delta_model'))=[];
conditions(strcmp(conditions,'delta_bin_midpoints'))=[];
mouse='15784';
mu={};
ymin={};
ymax={};

for i=1:length(conditions)
   mu{i}=mean(lag_mat.(conditions{i}).ref);
   ci=bootci(1e3,{@mean,lag_mat.(conditions{i}).ref},'type','per','alpha',.01);
   ymin{i}=ci(1,:)';
   ymax{i}=ci(2,:)';
   bootstr_mu{i}=mean(rnd_summary.(conditions{i}).ref);
   bootstr_min{i}=-3*std(rnd_summary.(conditions{i}).ref)';
   bootstr_max{i}=3*std(rnd_summary.(conditions{i}).ref)';
end

figs.changepoints=figure('PaperPositionMode','auto','Position',[100 100 800 250]);
whitebg(figs.changepoints);

for i=1:length(conditions)
    g(1,i)=gramm('x',([-90:90]/camera_fs)','y',{mu{i};bootstr_mu{i}},'ymin',{ymin{i};bootstr_min{i}},'ymax',{ymax{i};bootstr_max{i}},'color',1:2);
    g(1,i).axe_property('ylim',[-.05 .25],'xlim',[-3 3],'xtick',[-3 0 3],'ytick',[-.05:.05:.25],'TickDir','out','TickLength',[.015 .015]);
    g(1,i).geom_line;
    g(1,i).geom_interval('geom','area')
    g(1,i).no_legend;
    g(1,i).set_color_options('map',[1 0 0;.8 .8 .8]);
    g(1,i).set_title(conditions{i})
    g(1,i).set_names('x','','y','');
end

g.draw
% 
for i=1:length(g)
    offsetAxes(g(i).facet_axes_handles);
end
set(figs.changepoints,'color',[0 0 0]);
set(figs.changepoints,'color',[0 0 0],'InvertHardcopy','off')

% pdf is the only output that looks half decent here

markolab_multi_fig_save(figs.changepoints,'~/Desktop/quickfigs',[ mouse '_systemclub_changepoint_xcorr'] ,'png,fig,pdf','renderer','painters');

%%

% combine across mice?

files={};
files{1}='/Users/jmarkow/Desktop/workspace/photometry/15781/analysis/timescale_analysis.mat';
files{2}='/Users/jmarkow/Desktop/workspace/photometry/15783/analysis/timescale_analysis.mat';
files{3}='/Users/jmarkow/Desktop/workspace/photometry/15784/analysis/timescale_analysis.mat';

max_r_all={};
tmp1=[];

for i=1:length(files)
    load(files{i},'max_r','score_pks');
    max_r_all=[max_r_all max_r];
    tmp1=[tmp1;cellfun(@(x) median(diff(x)),score_pks)];
end

max_r_all=cat(1,max_r_all{:});
mins=min(max_r_all,[],2);
maxs=max(max_r_all,[],2);

max_r_all=(max_r_all-repmat(mins,[1 size(max_r_all,2)]))./(repmat(maxs,[1 size(max_r_all,2)])-repmat(mins,[1 size(max_r_all,2)]))
figs.timescale=figure('PaperPositionMode','auto','Position',[100 100 200 250]);
whitebg(figs.timescale);

tmp1=mean(tmp1);
tmp2=mean(max_r_all);
ci=bootci(1e3,{@mean,max_r_all},'type','per','alpha',.01);
g=gramm('x',tmp1*1/30,'y',tmp2,'ymin',ci(1,:),'ymax',ci(2,:));
g.geom_line;
g.geom_interval('geom','area');
g.no_legend;
g.set_names('x','','y','');
g.axe_property('ylim',[0 1],'ytick',[0 .5 1],'TickDir','out','TickLength',[.015 .015])
g.draw
offsetAxes(g.facet_axes_handles);

set(figs.timescale,'color',[0 0 0],'InvertHardcopy','off');


%%

% pca





