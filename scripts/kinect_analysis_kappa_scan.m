%%%
%%

% use two metrics, (1) average correlation between rps and model
% changepoints (2) difference in distributions

[nparams,nsessions]=size(labels);
peakthresh=.1;
peak_dist=cell(size(labels));

timer_upd=kinect_extract.proc_timer(nparams*nsessions,'frequency',10);
counter=0;
test=struct() 
tmp_corr=nan(nparams,nsessions);
kld=@(x,y) sum(x.*log2(x./y)); 
bins=[0:.005:2.5];

test.isi_rps=[];
for i=1:nsessions
    [~,locs]=findpeaks(zscore(extract_object(i).projections.rp_changepoint_score),'minpeakheight',.5);
    isi_rps=diff(locs)/30;
    test.isi_rps=[test.isi_rps;isi_rps];
end

for i=1:nparams
    
    test.isi{i}=[];
    
    for j=1:nsessions

        if ~isnan(labels{i,j})
            model_changepoints=abs(diff([-1;labels{i,j}(:)]))>0;
        else
            continue;
        end

        kernel=normpdf(-10:10,0,2);
        kernel=kernel./sum(kernel);
        
        tmp=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),...
            zscore(extract_object(j).projections.rp_changepoint_score));
        tmp_corr(i,j)=tmp(2,1);
        
        idxs=find(abs(diff([-1;labels{i,j}(:)]))>0);       
        isi=diff(idxs)/30;
        test.isi{i}=[test.isi{i};isi];
       
        counter=counter+1;
        timer_upd(counter);
        
    end
end


%%

kappas=cellfun(@(x) x.kappa,scan_dicts);
gammas=cellfun(@(x) x.gamma,scan_dicts);
nkappas=length(unique(kappas));
ngammas=length(unique(gammas));

measures={@mode,@mean,@std,@skewness};
measure_names={'mode','mean','std','skew'};

sum_stats=struct();

for i=1:length(measures)
    sum_stats.isi.(measure_names{i})=reshape(cellfun(measures{i},test.isi),ngammas,nkappas);
    sum_stats.isi_rps.(measure_names{i})=measures{i}(test.isi_rps);
end    

sum_stats.isi.corr=mean(reshape(tmp_corr,ngammas,nkappas,[]),3);

%%

sum_stats.isi.jsd=zeros(length(test.isi),1);

for i=1:length(test.isi)
   p1=histc(test.isi{i},bins);
   p2=histc(isi_rps,bins);
   
   if isempty(p1) | isempty(p2)
       continue;
   end
   
   p1=p1+eps;
   p2=p2+eps;
   p1=p1./sum(p1);
   p2=p2./sum(p2);
   m=(p1+p2)*.5;
  
   kld_use=~(p1==0|p2==0);
   
   p1=p1(kld_use);
   p2=p2(kld_use);
   sum_stats.isi.jsd(i)=.5*kld(p1,m)+.5*kld(p2,m);    
   
end

sum_stats.isi.jsd=reshape(sum_stats.isi.jsd,ngammas,nkappas);

sum_stats.isi.nstates=nan(nparams,1);

for i=1:nparams
    tmp=cat(2,labels{i,:});
    tmp(isnan(tmp))=[];
    sum_stats.isi.nstates(i)=length(unique(tmp));
end

sum_stats.isi.nstates=reshape(sum_stats.isi.nstates,ngammas,nkappas);

%% base on neural activity

phot=extract_object.get_photometry;
%[nkappas,nsessions]=size(labels);
nrnds=1e3;

gcamp_model_corr_rnd=zeros(nkappas,nsessions,nrnds);
rcamp_model_corr_rnd=zeros(nkappas,nsessions,nrnds);
gcamp_model_corr=zeros(nkappas,nsessions);
rcamp_model_corr=zeros(nkappas,nsessions);
cp_corr=zeros(nkappas,nsessions);
counter=0;

delta_win=10;
kernel=normpdf(-5:5,0,1);
kernel=kernel./sum(kernel);
rcamp_kernel=(ones(1,500)*.65).^[1:500];
rcamp_kernel=rcamp_kernel./sum(rcamp_kernel);
trace_use='dff_delta';

gcamp_score={};
rcamp_score={};
for i=1:nsessions
    
    % scores for correlation
    if isempty(phot(i).traces)
        continue;
    end
    
    gcamp_score{i}=markolab_deltacoef(zscore(phot(i).traces(1).dff),delta_win);
    smooth_rcamp=[conv(phot(i).traces(2).dff,rcamp_kernel(end:-1:1),'valid');zeros(numel(rcamp_kernel)-1,1)];
    rcamp_score{i}=markolab_deltacoef(zscore(smooth_rcamp),delta_win);
    
    gcamp_score{i}(gcamp_score{i}<=0|isnan(gcamp_score{i}))=0;
    rcamp_score{i}(rcamp_score{i}<=0|isnan(rcamp_score{i}))=0;
    
end

%%

% make a plot of how we get to the score fool

fig.scoresample=figure();

smooth_rcamp=[conv(phot(5).traces(2).dff,rcamp_kernel(end:-1:1),'valid');zeros(numel(rcamp_kernel)-1,1)];
ax=[];
ax(1)=subplot(2,2,1);
plot(phot(5).traces(1).dff);

ax(2)=subplot(2,2,2);
plot(smooth_rcamp);

ax(3)=subplot(2,2,3);
plot(gcamp_score{5});

ax(4)=subplot(2,2,4);
plot(rcamp_score{5});

linkaxes(ax,'x');

%%

timer_upd=kinect_extract.proc_timer(nparams*nsessions,'frequency',20);
counter=1

for i=1:nparams
    
    for j=1:nsessions
        
        counter=counter+1;
        
        if ~isnan(labels{i,j})
            model_changepoints=abs(diff([-1;extract_object(j).get_original_timebase(labels{i,j}(:))]))>0;
        else
            continue;
        end
        
        if isempty(phot(j).traces)
            continue;
        end
        
        %         tmp=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(1).(trace_use),phot_kernel,'same')));
        %         tmp2=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(2).(trace_use),phot_kernel,'same')));
        %
        
        wins=kinect_extract.window_data(zscore(gcamp_score{j}),find(model_changepoints)+1,5);
        wins2=kinect_extract.window_data(zscore(rcamp_score{j}),find(model_changepoints)+1,5);
        
        gcamp_tstat=mean(nanstd(wins,[],2))./sqrt(size(wins,2)); 
        rcamp_tstat=mean(nanstd(wins2,[],2))./sqrt(size(wins2,2));
        
        gcamp_model_corr(i,j)=max(nanmean(wins,2))./gcamp_tstat;
        rcamp_model_corr(i,j)=max(nanmean(wins2,2))./rcamp_tstat;
        
        %         wins3=kinect_extract.window_data(zscore(extract_object(j).projections.rp_changepoint_score),find(model_changepoints)+1,5);
        %         cp_corr(i,j)=max(nanmean(wins3,2));
        
        %
        
%         gcamp_model_corr(i,j)=tmp(2,1);
%         rcamp_model_corr(i,j)=tmp2(2,1);
        
        %
        timer_upd(counter);
        
        
    end
    
end

gcamp_model_corr=reshape(gcamp_model_corr,nkappas,ngammas,[]);
rcamp_model_corr=reshape(rcamp_model_corr,nkappas,ngammas,[]);

%%

xvec=log10(unique(kappas));
yvec=log10(unique(gammas));
plot_stats={'mean','mode','cv','std','skew','nstates','jsd','corr'};
figs.summary_stats=figure();
sum_stats.isi.cv=sum_stats.isi.std./sum_stats.isi.mean;
[ynans,xnans]=find(isnan(sum_stats.isi.mean));

for i=1:length(plot_stats)
   subplot(3,3,i);
   
   clims(1)=min(sum_stats.isi.(plot_stats{i})(:));
   clims(2)=max(sum_stats.isi.(plot_stats{i})(:))
   new_clims=clims*10;
   new_clims(1)=floor(new_clims(1))/10;
   new_clims(2)=ceil(new_clims(2))/10;
   
  
   plotim=sum_stats.isi.(plot_stats{i});
   plotim=(plotim-new_clims(1))./(range(new_clims));
   plotim(plotim<0)=0;
   plotim(plotim>1)=1;
   
   rgbim=ind2rgb(uint16(plotim*size(fee_map,1)),fee_map);
   
   for j=1:length(xnans)
       rgbim(ynans(j),xnans(j),:)=[.7 .7 .7];
   end
    
   image(xvec,yvec,rgbim);
   axis xy;
   
   c=colorbar();
   colormap(fee_map)
   yticks=get(c,'ticks');
   box off;
   set(gca,'TickDir','out','TickLength',[.015 .015],'XTick',[0:5:xvec(end)],'YTick',[0:5:yvec(end)]);
    set(c,'ticks',[0 1],'ticklabels',[new_clims]);
   title(plot_stats{i});
end

%%

fig.d1d2=figure();

subplot(1,2,1);
   
useim=mean(gcamp_model_corr(:,:,3:6),3);
clims(1)=min(useim(:));
clims(2)=max(useim(:));
new_clims=clims*10;
new_clims(1)=floor(new_clims(1))/10;
new_clims(2)=ceil(new_clims(2))/10;

plotim=useim;
plotim=(plotim-new_clims(1))./(range(new_clims));
plotim(plotim<0)=0;
plotim(plotim>1)=1;

rgbim=ind2rgb(uint16(plotim*size(fee_map,1)),fee_map);

for j=1:length(xnans)
   rgbim(ynans(j),xnans(j),:)=[.7 .7 .7];
end

image(xvec,yvec,rgbim);
axis xy;
 
c=colorbar();
colormap(fee_map)
yticks=get(c,'ticks');
box off;
set(gca,'TickDir','out','TickLength',[.015 .015],'XTick',[0:5:xvec(end)],'YTick',[0:5:yvec(end)]);
set(c,'ticks',[0 1],'ticklabels',[new_clims]);
title('D2');


subplot(1,2,2);
   
useim=mean(rcamp_model_corr(:,:,3:6),3);
clims(1)=min(useim(:));
clims(2)=max(useim(:));
new_clims=clims*10;
new_clims(1)=floor(new_clims(1))/10;
new_clims(2)=ceil(new_clims(2))/10;

plotim=useim;
plotim=(plotim-new_clims(1))./(range(new_clims));
plotim(plotim<0)=0;
plotim(plotim>1)=1;

rgbim=ind2rgb(uint16(plotim*size(fee_map,1)),fee_map);

for j=1:length(xnans)
   rgbim(ynans(j),xnans(j),:)=[.7 .7 .7];
end

image(xvec,yvec,rgbim);
axis xy;

c=colorbar();
colormap(fee_map)
yticks=get(c,'ticks');
box off;
set(gca,'TickDir','out','TickLength',[.015 .015],'XTick',[0:5:xvec(end)],'YTick',[0:5:yvec(end)]);
set(c,'ticks',[0 1],'ticklabels',[new_clims]);
title('D1');

% TODO:  add the gradient here, just smooth and look at fx, also need to
% pick right timescale for deltascores

%%



%%



% %%
% counter=0;
% 
% timer_upd=kinect_extract.proc_timer(nkappas*nrnds,'frequency',20);
% 
% for ii=1:nrnds
%     
%     for i=1:nkappas
%         
%         for j=1:nsessions
%             
%             if ~isnan(labels{i,j})
%                 model_changepoints=abs(diff([-1;extract_object(j).get_original_timebase(labels{i,j})]))>0;
%             else
%                 disp('test');
%                 continue;
%             end
% 
%             if isempty(phot(j).traces)
%                 continue;
%             end
%                        
%             tmp=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(1).dff_rnd(:,ii),phot_kernel,'same')));
%             tmp2=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(2).dff_rnd(:,ii),phot_kernel,'same')));
% 
%             gcamp_model_corr_rnd(i,j,ii)=tmp(2,1);
%             rcamp_model_corr_rnd(i,j,ii)=tmp2(2,1);
%                       
%         end
%         
%         counter=counter+1;
%         timer_upd(counter);
%         
%     end
%     
% end
% 
% save('_analysis/kappa_gamma_analysis.mat','gcamp_model_corr','rcamp_model_corr','gcamp_model_corr_rnd','rcamp_model_corr_rnd','-v7.3');

% %%
% 
% % plot the em effin thing
% 
% edge=1:size(test,1);
% figure();
% imagesc(edge,edge,mean(imgaussfilt(test(:,:,4:6),.5),3))
% %imagesc(edge,edge,mean(test(:,:,6),3))
% colormap(jet)
% c=colorbar()
% ylabel(c,'Model/Photometry correlation','FontSize',12)
% xlabel('Kappa')
% ylabel('Gamma')
% set(gca,'YTick',[2:2:14],'YTickLabel',{'10^2','10^4','10^6','10^8','10^{10}','10^{12}','10^{14}'})
% set(gca,'XTick',[2:2:14],'XTickLabel',{'10^2','10^4','10^6','10^8','10^{10}','10^{12}','10^{14}'})
% set(gca,'TickLength',[.015 .015])
% 
% set(gca,'TickDir','out')
% box off
% 
% %%
% 
% use_session=3:6;
% 
% use_rnd=rcamp_model_corr_rnd(:,use_session,:);
% rcamp_model_corr_z=(mean(rcamp_model_corr(:,3:6)')-mean(squeeze(mean(use_rnd,2))'))./std(squeeze(mean(use_rnd,2))');
% 
% use_rnd=gcamp_model_corr_rnd(:,use_session,:);
% gcamp_model_corr_z=(mean(gcamp_model_corr(:,3:6)')-mean(squeeze(mean(use_rnd,2))'))./std(squeeze(mean(use_rnd,2))');


