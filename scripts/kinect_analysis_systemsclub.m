%
% assumes masked data and bounded data are loaded in;

% generate a simple movie that rotates to demonstrate 
% nature of data


[height,width,nframes]=size(depth_bounded_rotated);
frames=single(reshape(depth_bounded_rotated.*int16(depth_bounded_cable_mask_rotated),height*width,[]));
rps=kinect_gaussproj(frames',600);
rps=zscore(zscore(rps)');
disp_rps=imgaussfilt(filter(ones(5,1)/5,1,rps')',1.5);
scores=frames'*features.svd.U;

%%

% animate that shit

v=VideoWriter('kinect_example.mp4','mpeg-4');
v.FrameRate=30;
v.Quality=100;

frames=5e3:10.1e3;

% how long to rotate for

rotate_frames=300;
start=300;
rps_width=300;
az=ones(1,length(frames))*0;
el=ones(1,length(frames))*90;
az(start+1:start+rotate_frames)=linspace(0,-15,rotate_frames);
el(start+1:start+rotate_frames)=linspace(90,13,rotate_frames);
az(start+rotate_frames+1:start+rotate_frames*2)=linspace(-15,0,rotate_frames);
el(start+rotate_frames+1:start+rotate_frames*2)=linspace(13,90,rotate_frames);

fig=figure();
whitebg(fig,[0 0 0]);
set(fig,'color',[0 0 0]);
set(fig,'position',[100 100 1e3 500 ]);
[xx,yy]=meshgrid(1:512,1:424);
zz=rand(size(depth_masked(:,:,1)));
ax=axes('ydir','rev','units','pixels','position',[50 50 512/1 424/1]);
h=surf(xx,yy,zz,'parent',ax);shading interp
hold on;
axis(ax,[0 512 0 424 0 80],'off');caxis([0 40]);
l1=plot3([100 100],[424 424],[0 10],'w-','parent',ax,'linewidth',1.5);
l2=plot3([100 140],[424 424],[0 0],'w-','parent',ax,'linewidth',1.5);
l3=plot3([100 100],[384 424],[0 0],'w-','parent',ax,'linewidth',1.5);
colormap(bone);

ax2=axes('ydir','rev','units','pixels','position',[650 474-160 160 160]);
h2=imagesc(rand(size(depth_bounded_rotated(:,:,1))),'parent',ax2);caxis([0 40]);
axis(ax2,'off');

ax3=axes('units','pixels','position',[520 200 400 100]);
h3=imagesc(disp_rps(:,frames(1):frames(1)+rps_width));
h4=patch([0 rps_width+1 rps_width+1 0],[ 1 1 size(rps,1) size(rps,1) ],0,'facecolor',[0 0 0],'edgecolor','none');
caxis([-1 1]);
axis(ax3,'off');
set(ax3,'xlim',[1 rps_width]);

score_plot=filter(ones(5,1)/5,1,zscore(scores(:,:)))+repmat([size(scores,2):-1:1]*2,[size(scores,1) 1]);

ax4=axes('units','pixels','position',[520 70 400 100]);
h5=plot(score_plot(frames(1):frames(1)+rps_width,1:10),'linewidth',1.5);
ylimits=ylim(ax4);
h6=patch([0 rps_width+1 rps_width+1 0],[ ylimits(1) ylimits(1) ylimits(2) ylimits(2) ],0,'facecolor',[0 0 0],'edgecolor','none');
axis(ax4,'off');
set(ax4,'xlim',[1 rps_width]);

open(v);
timer_upd=kinect_proctimer(length(frames));


for i=frames
    set(h,'XData',xx,'YData',yy,'ZData',depth_masked(:,:,i));
    set(h2,'CData',depth_bounded_rotated(:,:,i));
    if (i-(frames(1)-1))<=rps_width    
        set(h4,'xdata',[i-(frames(1)-1) rps_width+1 rps_width+1 i-(frames(1)-1)]');
        set(h6,'xdata',[i-(frames(1)-1) rps_width+1 rps_width+1 i-(frames(1)-1)]');
    else
        set(h3,'cdata',disp_rps(:,i-rps_width:i));
        for j=1:length(h5)
            h5(j).YData=score_plot(i-rps_width:i,j);
        end
    end
    view(ax,az(i-(frames(1)-1)),el(i-(frames(1)-1)));
    im=getframe(fig);
    writeVideo(v,im.cdata);
    pause(eps);
    timer_upd(i-(frames(1)-1));
end

close(v);





 