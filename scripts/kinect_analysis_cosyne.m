
%%
% 
use_session=6;
phot=extract_object.get_photometry;
cp_score=extract_object(use_session).get_original_timebase(extract_object(use_session).projections.rp_changepoint_score);

[r1,lags]=xcorr(zscore(phot(use_session).traces(1).dff),zscore(cp_score),200,'coeff');
[r2,~]=xcorr(zscore(phot(use_session).traces(2).dff),zscore(cp_score),200,'coeff');

% wanna window the data bro?

gcamp_wins=kinect_extract.vec_to_mat(zscore(phot(use_session).traces(1).dff),5e2,250);
rcamp_wins=kinect_extract.vec_to_mat(zscore(phot(use_session).traces(2).dff),5e2,250);
cp_wins=kinect_extract.vec_to_mat(zscore(cp_score),5e2,250);
nwins=size(gcamp_wins,2);

max_lag=100;
r1_mat=zeros(max_lag*2+1,nwins);
r2_mat=zeros(max_lag*2+1,nwins);

for i=1:size(gcamp_wins,2)
    r1_mat(:,i)=xcorr(gcamp_wins(:,i),cp_wins(:,i),max_lag,'coeff');
    r2_mat(:,i)=xcorr(rcamp_wins(:,i),cp_wins(:,i),max_lag,'coeff');
end

r1_ci=bootci(1e3,{@mean,r1_mat'},'type','per');
r2_ci=bootci(1e3,{@mean,r2_mat'},'type','per');

r1_mu=mean(r1_mat,2);
r2_mu=mean(r2_mat,2);

r1_plot=[r1_ci(1,:)' r1_mu r1_ci(2,:)'];
r2_plot=[r2_ci(1,:)' r2_mu r2_ci(2,:)'];

% scale for display


%%

labels2=reshape(labels,sqrt(size(labels,1)),sqrt(size(labels,1)),[]);


cp_score=abs(diff([-1;extract_object(use_session).get_original_timebase(labels2{3,6,use_session})]))>0;
kernel=normpdf(-20:20,0,2);
cp_score=conv(single(cp_score),kernel,'same');

[r1,lags]=xcorr(zscore(phot(use_session).traces(1).dff),zscore(cp_score),200,'coeff');
[r2,~]=xcorr(zscore(phot(use_session).traces(2).dff),zscore(cp_score),200,'coeff');

%%

r1_plot=bsxfun(@minus,r1_plot,min(r1_plot(:)));
r1_plot=bsxfun(@rdivide,r1_plot,(max(r1_plot(:))-min(r1_plot(:))));

r2_plot=bsxfun(@minus,r2_plot,min(r2_plot(:)));
r2_plot=bsxfun(@rdivide,r2_plot,(max(r2_plot(:))-min(r2_plot(:))));

fig=figure();
markolab_shadeplot([-max_lag:max_lag]/30,[r2_plot(:,1) r2_plot(:,3)]','g','k',2);
hold on
plot([-max_lag:max_lag]/30,r2_plot(:,2),'k-');
markolab_shadeplot([-max_lag:max_lag]/30,[r1_plot(:,1) r1_plot(:,3)]','r','k',2);
hold on
plot([-max_lag:max_lag]/30,r1_plot(:,2),'k-');


%%

movie_session=6;
movie_frames=6e3:10e3;

nparams=size(labels,1);

labels2=reshape(labels,sqrt(nparams),sqrt(nparams),size(labels,2));
scan_dicts2=reshape(scan_dicts,sqrt(nparams),[]);

for i=1:sqrt(nparams)
   kappa=scan_dicts2{2,i}.kappa;
   extract_object.load_model_labels(labels2(2,i,:),[]);
   extract_object(6).animate_data(movie_frames,[],['test_movie_session6_30pcs_kappa1e' num2str(log10(kappa)) 'gamma1e3_frames6e3-10e3.mp4']);
end

%%

movie_session=6;
movie_frames=6e3:10e3;

nparams=size(labels,1);

labels2=reshape(labels,sqrt(nparams),sqrt(nparams),size(labels,2));
scan_dicts2=reshape(scan_dicts,sqrt(nparams),[]);

for i=1:sqrt(nparams)
   kappa=scan_dicts2{2,i}.kappa;
   extract_object.load_model_labels(labels2(3,i,:),[]);
   extract_object(6).animate_data(movie_frames,[],['test_movie_session6_kappa1e' num2str(log10(kappa)) 'gamma1e3_frames6e3-10e3.mp4']);
end