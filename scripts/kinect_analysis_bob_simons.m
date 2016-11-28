%%%

% take some shit, do some shit


% assume pcs are loaded in, use for reconstruction


% use day

use_session=7;

% reconstruct the frames?

load(metadata(use_session).frames_file,'depth_bounded_rotated','depth_bounded_cable_mask_rotated');
orientation=scalars(use_session).orientation;
centroid=scalars(use_session).centroid;
%%
% threshold data

frame_data=single(depth_bounded_rotated);
frame_data(log(depth_bounded_cable_mask_rotated)<-11)=nan;
frame_data(depth_bounded_rotated<10)=0;
frame_data=reshape(frame_data,80^2,[]);

frame_data=frame_data';
replace_idx=isnan(frame_data);

for i=1:size(frame_data,1)
  idx=isnan(frame_data(i,:));
  this_frame=frame_data(i,:);
  frame_data(i,idx)=mean(this_frame(this_frame>0));
end

iters=20;

for i=1:iters
  i
  mu=mean(frame_data,1);
  scores=bsxfun(@minus,frame_data,mu)*features.pca.coeffs(:,1:50);
  recon=scores(:,1:20)*features.pca.coeffs(:,1:20)';
  recon=bsxfun(@plus,recon,mu);
  frame_data(replace_idx)=recon(replace_idx);
end

frame_data=reshape(frame_data',80,80,[]);
recon=reshape(recon',80,80,[]);

%%

% use the reconstruciton to do stuff

spines=imagesc(squeeze(mean(frame_data(40:60,:))));




%%
% now either put recon back in

full_recon=zeros(424,512,size(frame_data,3),'int16');
box_size=[80 80];

for i=1:size(frame_data,3)

  insert_mouse=recon(:,:,i);
  insert_mouse=imgaussfilt(imrotate(insert_mouse,orientation(i),'bilinear','crop'),.1);


  coords_x=centroid(i,1)-(box_size(2)/2-1):centroid(i,1)+box_size(2)/2;
  coords_y=centroid(i,2)-(box_size(1)/2-1):centroid(i,2)+box_size(1)/2;

  full_recon(coords_y,coords_x,i)=insert_mouse;
  figure(1);
  imagesc(full_recon(:,:,i));caxis([10 80]);
  pause(.01);

end


%%
