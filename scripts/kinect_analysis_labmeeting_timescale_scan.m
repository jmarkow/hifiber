%% cycle over different windows sizes

win_size=[2:4 5:5:40];
bins=[0:100:2e3];
bin_density={};
delta_score={};
obs_r={};
ibi={};
bootvals={};

% get delta coefficients

for i=1:length(win_size)
    
    fprintf(1,'%g\n',win_size(i));
    rps_deltas=markolab_deltacoef(rps,win_size(i)); % delta coefficients, lag of 4 frames
    delta_score{i}=sum(abs(rps_deltas)>.1); % binarize deltas
    h=normpdf([-10:10],0,.5); % gauss smooth
    delta_score{i}=conv(delta_score{i},h,'same'); % convolve

    thresh=mean(delta_score{i})+.5*std(delta_score{i});
    idx=1:length(delta_score{i})-1;

    %block_onsets=delta_score{i}(idx)<thresh&delta_score{i}(idx+1)>thresh;
    [block_peaks,block_locs]=findpeaks(delta_score{i},'minpeakheight',thresh);
    ibi{i}=diff((block_locs)*ms_per_frame);
    
    nboots=1e3;
    bootvals_max=nan(1,nboots);
    bootvals_min=nan(1,nboots);

    for j=1:nboots
        scr=markolab_phase_scramble_1d(zscore(gcamp_interp),0);
        r=xcorr(zscore(scr),zscore(delta_score{i}),25,'coeff');
        bootvals_max(j)=max(r);
        bootvals_min(j)=min(r);
    end
    
    bootvals{i}=bootvals_max;
    bin_density{i}=histc(ibi{i},bins);
    [obs_r{i},lags]=xcorr(zscore(gcamp_interp),zscore(delta_score{i}),25,'coeff');

end

%%

figure();

bootmu=cellfun(@mean,bootvals);
bootstd=cellfun(@mean,bootvals);
corrvals=cellfun(@max,obs_r);
med_ibi=cellfun(@median,ibi);


subplot(1,2,1);
yyaxis left;


plot((corrvals-bootmu)./bootstd,'b.-','markersize',15);
ylabel('Max r (Z, bootstrap)');
yyaxis right;
plot(med_ibi,'r.-','markersize',15);
grid minor;
set(gca,'MinorGridLineStyle','-');
box off;


subplot(1,2,2);

yyaxis left;
plot((corrvals-bootmu)./bootstd,'b.-','markersize',15);

yyaxis right;
plot(med_ibi,'r.-','markersize',15);
grid minor;
set(gca,'MinorGridLineStyle','-');
box off;
ylabel('Med. block duration (ms)');

ylim([0 1e3]);
