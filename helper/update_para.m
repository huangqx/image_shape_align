function [Para] = update_para(Para)

for i=1:1600
    h=figure('visible','off');
    disp(i);
    set(h,'Position',[1 1 i 10])   
    f=getframe;
    [image, map] = frame2im(f);
    Para.rWidths(1,i)=size(image,2);
    close(h)
end

for i=1:1600
    h=figure('visible','off');
    disp(i);
    set(h,'Position',[1 1 10 i])   
    f=getframe;
    [image, map] = frame2im(f);
    Para.rHeights(1,i)=size(image,1);
    close(h)
end


end