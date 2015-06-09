function res = track_video2(frames,m_average,x_average,sizes,movieName,frameName)

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
for i=1:1:size(m_average,3)
    hold off
    imshow(frames(:,:,i),[lb ub],'Border','tight')
    hold on
    axis off
    m=m+1;
    % set(gca, 'LooseInset', get(gca,'TightInset'))

   plot([m_average(1:sizes,1,i);m_average(1,1,i)],[m_average(1:sizes,2,i);m_average(1,2,i)],'Linewidth',2,'Color',[0.99 0.99 0.99])
   plot([m_average(sizes+1:end,1,i);m_average(sizes+1,1,i)],[m_average(sizes+1:end,2,i);m_average(sizes+1,2,i)],'Linewidth',2,'Color',[0.99 0.99 0.99])
   if (~isempty(x_average))
    plot(x_average(:,1,i),x_average(:,2,i),'o','Linewidth',2,'Color',[0.99 0.99 0.99])
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
