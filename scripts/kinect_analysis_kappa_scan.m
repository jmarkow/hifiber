%%%
%%
%
%

% use two metrics, (1) average correlation between rps and model
% changepoints (2) difference in distributions

[nkappas,nsessions]=size(labels);
peakthresh=.1;
peak_dist=cell(size(labels));

timer_upd=kinect_extract.proc_timer(nkappas*nsessions,'frequency',10);
counter=0;
test=struct() 

kld=@(x,y) sum(x.*log2(x./y)); 
bins=[0:.01:5];

for i=1:nkappas
    
    
    for j=1:nsessions
        
        %        [~,locs]=findpeaks(extract_object(j).projections.rp_changepoint_score,'minpeakheight',peakthresh);
        %
        if ~isnan(labels{i,j})
            model_changepoints=abs(diff([-1;labels{i,j}(:)]))>0;
        else
            continue;
        end
        %
        %        peak_dist{i,j}=zeros(size(model_changepoints));
        %
        %        for k=1:length(model_changepoints)
        %           peak_dist{i,j}(k)=min(abs(model_changepoints(k)-locs));
        %        end
        
        kernel=normpdf(-10:10,0,2);
        kernel=kernel./sum(kernel);
        
        tmp=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(extract_object(j).projections.rp_changepoint_score));
        
        test.corr(i,j)=tmp(2,1);
        
        idxs=find(abs(diff([-1;labels{i,j}(:)]))>0);
        [~,locs]=findpeaks(zscore(extract_object(j).projections.rp_changepoint_score),'minpeakheight',.25);
        
        isi=diff(idxs)/30;
        isi_rps=diff(locs)/30;
        
        p1=histc(isi,bins);
        p2=histc(isi_rps,bins);
        
        p1=p1+eps;
        p2=p2+eps;
        p1=p1./sum(p1);
        p2=p2./sum(p2);
        m=(p1+p2)*.5;
        
        %kld_use=~(p1==0|p2==0);
        
        %p1=p1(kld_use);
        %p2=p2(kld_use);
        
        test.mu(i,j)=mean(isi);
        test.std(i,j)=std(isi);
        test.skew(i,j)=skewness(isi);
        test.mode(i,j)=mode(isi);
        
        test.jsd(i,j)=.5*kld(p1,m)+.5*kld(p2,m);
        
        counter=counter+1;
        timer_upd(counter);
        
    end
end

kappas=cellfun(@(x) x.kappa,scan_dicts);
gammas=cellfun(@(x) x.gamma,scan_dicts);

quantities=fieldnames(test);
edge_size=sqrt(nkappas);

for i=1:length(quantities)
    test.(quantities{i})=reshape(test.(quantities{i}),length(unique(gammas)),length(unique(kappas)),[]);
end

%%

% make a series of panels fool

%% base on neural activity

phot=extract_object.get_photometry;
[nkappas,nsessions]=size(labels);
nrnds=1e3;
 
gcamp_model_corr_rnd=zeros(nkappas,nsessions,nrnds);
rcamp_model_corr_rnd=zeros(nkappas,nsessions,nrnds);
gcamp_model_corr=zeros(nkappas,nsessions);
rcamp_model_corr=zeros(nkappas,nsessions);
cp_corr=zeros(nkappas,nsessions);
counter=0;

kernel=normpdf(-5:5,0,1);
kernel=kernel./sum(kernel);
phot_kernel=(ones(1,10)*.5).^[1:10];
phot_kernel=phot_kernel./sum(phot_kernel);
timer_upd=kinect_extract.proc_timer(nkappas*nsessions,'frequency',20);

for i=1:nkappas
    
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
        
                tmp=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(1).dff,phot_kernel,'same')));
                tmp2=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(2).dff,phot_kernel,'same')));
        
        
%         wins=kinect_extract.window_data(zscore(phot(j).traces(1).dff),find(model_changepoints)+1,5);
%         wins2=kinect_extract.window_data(zscore(phot(j).traces(2).dff),find(model_changepoints)+1,5);
%         
%         gcamp_model_corr(i,j)=max(nanmean(wins,2));
%         rcamp_model_corr(i,j)=max(nanmean(wins2,2));
%         wins3=kinect_extract.window_data(zscore(extract_object(j).projections.rp_changepoint_score),find(model_changepoints)+1,5);
%         cp_corr(i,j)=max(nanmean(wins3,2));
        
%         
        gcamp_model_corr(i,j)=tmp(2,1);
        rcamp_model_corr(i,j)=tmp2(2,1);
%         
        timer_upd(counter);
        
        
    end
    
end


%%
counter=0;

timer_upd=kinect_extract.proc_timer(nkappas*nrnds,'frequency',20);

for ii=1:nrnds
    
    for i=1:nkappas
        
        for j=1:nsessions
            
            if ~isnan(labels{i,j})
                model_changepoints=abs(diff([-1;extract_object(j).get_original_timebase(labels{i,j})]))>0;
            else
                disp('test');
                continue;
            end

            if isempty(phot(j).traces)
                continue;
            end
                       
            tmp=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(1).dff_rnd(:,ii),phot_kernel,'same')));
            tmp2=corrcoef(zscore(conv(single(model_changepoints),kernel,'same')),zscore(conv(phot(j).traces(2).dff_rnd(:,ii),phot_kernel,'same')));

            gcamp_model_corr_rnd(i,j,ii)=tmp(2,1);
            rcamp_model_corr_rnd(i,j,ii)=tmp2(2,1);
                      
        end
        
        counter=counter+1;
        timer_upd(counter);
        
    end
    
end

save('_analysis/kappa_gamma_analysis.mat','gcamp_model_corr','rcamp_model_corr','gcamp_model_corr_rnd','rcamp_model_corr_rnd','-v7.3');

%%

% plot the em effin thing

edge=1:size(test,1);
figure();
imagesc(edge,edge,mean(imgaussfilt(test(:,:,4:6),.5),3))
%imagesc(edge,edge,mean(test(:,:,6),3))
colormap(jet)
c=colorbar()
ylabel(c,'Model/Photometry correlation','FontSize',12)
xlabel('Kappa')
ylabel('Gamma')
set(gca,'YTick',[2:2:14],'YTickLabel',{'10^2','10^4','10^6','10^8','10^{10}','10^{12}','10^{14}'})
set(gca,'XTick',[2:2:14],'XTickLabel',{'10^2','10^4','10^6','10^8','10^{10}','10^{12}','10^{14}'})
set(gca,'TickLength',[.015 .015])

set(gca,'TickDir','out')
box off

%%

use_session=3:6;

use_rnd=rcamp_model_corr_rnd(:,use_session,:);
rcamp_model_corr_z=(mean(rcamp_model_corr(:,3:6)')-mean(squeeze(mean(use_rnd,2))'))./std(squeeze(mean(use_rnd,2))');

use_rnd=gcamp_model_corr_rnd(:,use_session,:);
gcamp_model_corr_z=(mean(gcamp_model_corr(:,3:6)')-mean(squeeze(mean(use_rnd,2))'))./std(squeeze(mean(use_rnd,2))');


