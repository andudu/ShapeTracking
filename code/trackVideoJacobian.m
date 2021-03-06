 
function res = trackVideoJacobian(frames,m_average,x_average,alpha_average,sizes,movieName,frameName)

% function res = track_video(frames,m_average,x_average,movieName,frameName)


if (~isempty(movieName))
  mov = avifile([movieName,'.avi'], 'COMPRESSION', 'None', 'fps', 5,'VideoName', movieName);
end
lb = min(min(min(frames)));
ub = max(max(max(frames)));

figure(2);
colormap('gray')
cla

hold off
j=0;
m=1;

[pts_X, pts_Y] = meshgrid(1:5:size(frames,2),1:5:size(frames,1));
pts_new = [pts_X(:) pts_Y(:)];

x_new = x_average(:,:,1);

size(m_average)
size(frames)
sizes(1,1)

L_out = roipoly(frames(:,:,1),m_average(1:sizes(1),1,1),m_average(1:sizes(1),2,1));
L_in = roipoly(frames(:,:,1),m_average(sizes(1)+1:end,1,1),m_average(sizes(1)+1:end,2,1));
tube = L_out.*(1-L_in);



ind = find(tube>0); 
[I_idx,J_idx] = ind2sub(size(tube),ind);
pts = [J_idx I_idx];

pts_new = pts;
J = ones(size(pts_new,1),1);

for i=1:1:size(m_average,3)-1

    % ss[pts_new x_new,J_new] = updateJacobian([pts_X(:) pts_Y(:)],x_average(:,:,i),30,0.1,10,alpha_average(:,:,i));
    J_old = J;
    [pts_new x_new,J] = updateJacobian(pts_new,x_average(:,:,i),30,0.1,10,alpha_average(:,:,i));

    J = J.*J_old;
    
    max(J)
    min(J)
    hold off
    % imagesc(reshape(J_new,size(frames,1),size(frames,2)))
    %contour(pts_new(:,1),pts_new(:,2),J_new)
    % plot3(pts_new(:,1),pts_new(:,2),J_new,'r.')
    [X,Y] = meshgrid(1:1:size(frames,2),1:1:size(frames,1));
    Z = griddata(pts_new(:,1),pts_new(:,2),J,X,Y,'nearest');

    L_out = roipoly(frames(:,:,1),m_average(1:sizes(1),1,i+1),m_average(1:sizes(1),2,i+1));
    L_in = roipoly(frames(:,:,1),m_average(sizes(1)+1:end,1,i+1),m_average(sizes(1)+1:end,2,i+1));
    tube = L_out.*(1-L_in);

    imshow(Z,[0 2],'Border','tight')
    pos = get(gca,'Position');

    if (i==i)
      colorbar('location','West','Position',[pos(1)+0.05  pos(2) .02 pos(4)])
      % cpos =  get(colorbar,'Position')
      % cpos(3) = cpos(3)/4;
      % colorbar('off')
      % colorbar('Position',cpos,'location','WestOutside')
      % set(colorbar,'Position',cpos)
    else
      % colorbar('location','WestOutside','Position',cpos)
    end
    set(gca,'Position',pos)
    %% cpos(4) = cpos(4)/2;
    %set(colorbar,'Position',cpos)
 
    colormap(jet)
    hold on ;
    plot(pts_new(:,1),pts_new(:,2),'k.')
    %plot(x_new(:,1),x_new(:,2),'bo')
    %plot(x_average(:,1,i+1),x_average(:,2,i+1),'go')
    %plot(m_average(:,1,i+1),m_average(:,2,i+1),'g.','Linewidth',2)

    %imshow(frames(:,:,i),[lb ub],'Border','tight')
    hold on
    axis off
    m=m+1;
    % set(gca, 'LooseInset', get(gca,'TightInset'))

   plot([m_average(1:sizes,1,i+1);m_average(1,1,i+1)],[m_average(1:sizes,2,i+1);m_average(1,2,i+1)],'Linewidth',2,'Color',[0.99 0.99 0.99])
   plot([m_average(sizes+1:end,1,i+1);m_average(sizes+1,1,i+1)],[m_average(sizes+1:end,2,i+1);m_average(sizes+1,2,i+1)],'Linewidth',2,'Color',[0.99 0.99 0.99])
   if (~isempty(x_average))
    plot(x_average(:,1,i+1),x_average(:,2,i+1),'o','Linewidth',2,'Color',[0.99 0.99 0.99])
   end
   axis equal
   pause(0.05)

    if(~isempty(movieName))
       % frm = getframe(gcf);     mov = addframe(mov,frm) ;
       j=j+1;
       name = sprintf('%s%d',frameName,j);
       print ('-dpng', name)
       print ('-dpsc2', name)
       saveas (gcf,name,'fig')
       im = imread([name '.png']);
       frm = im2frame(im);
       mov = addframe(mov,frm);
       if (isempty(frameName))
	  delete([name '.png']);
          delete([name '.fig']);
          delete([name '.ps']);
       end       
    end
    if (~isempty(frameName)&&isempty(movieName))
      j=j+1;
      name = sprintf('%s%d',frameName,j);
      print ('-dpng' ,name)
      print ('-dpsc2', name)
      saveas (gcf,name,'fig')
    end


end


if (~isempty(movieName))
  mov = close(mov);
end
res = 0;
