 
function res = track_video_MCMC2(frames,m_average,x_average,sizes,movieName,frameName)

% res = track_video_MCMC(frames,m_average,x_average,movieName,frameName)



if (~isempty(movieName))
  mov = avifile([movieName,'.avi'], 'COMPRESSION', 'None', 'fps', 2,'VideoName', movieName);
end

lb = min(min(min(frames)));
ub = max(max(max(frames)));

% making sure there are no holes in the curve
m_average = [cat(1,m_average(1:sizes,:,:,:),m_average(1,:,:,:)); cat(1,m_average(sizes+1:end,:,:,:),m_average(sizes+1,:,:,:))];


% f = figure('visible','off');
f = figure(1);
colormap('gray')
cla

hold off
j=0;
m=1;

for i=1:1:size(m_average,4)
    hold off
    imshow(frames(:,:,i),[lb ub],'Border','tight')
    %imagesc(frames(:,:,i),[lb ub])
    hold on
    axis off
    m=m+1;

 plot(squeeze(m_average(1:sizes+1,1,:,i)),squeeze(m_average(1:sizes+1,2,:,i)),'Linewidth',2)
 plot(squeeze(m_average(sizes+2:end,1,:,i)),squeeze(m_average(sizes+2:end,2,:,i)),'Linewidth',2)
 plot(squeeze(x_average(:,1,:,i)),squeeze(x_average(:,2,:,i)),'o','Linewidth',2)
% axis equal 
 pause(0.01)
    
    
    if(~isempty(movieName))
       % frm = getframe(gcf);     mov = addframe(mov,frm) ;
       j=j+1;
       name = sprintf('%s%d',frameName,j);
       print ('-dpng',name)
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
      saveas(gcf,name,'png')
      saveas (gcf,name,'fig')
      print ('-dpsc2', name)
    end
end


if (~isempty(movieName))
  mov = close(mov);
end
res = 0;
set(f,'visible','on')
