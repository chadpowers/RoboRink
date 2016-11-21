function [out] = hsvRange(file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

swingObj = VideoReader(file);
vidWidth = swingObj.Width;
vidHeight = swingObj.Height;
scale = 1/5;
global a b c next shift

a = [0,1];
b = [0,1];
c = [0,1];
shift = 1;

%Create a movie structure array, mov.
mov = struct('cdata',zeros(vidHeight*scale,vidWidth*scale,3,'uint8'),...
    'colormap',[]);

hf = figure(3);
set(hf,'position',[10 10 vidWidth/2 vidHeight/2]);
set (hf,  'KeyPressFcn', @(src,eventdata)keyDownListener(src,eventdata));
set (hf,'KeyReleaseFcn', @(src,eventdata)keyUpListener(src,eventdata));

k = 1;
while hasFrame(swingObj) && ishandle(hf)
    values = readFrame(swingObj);
    d = rgb2hsv(values);
    d = imresize(d,scale);
    next = 0;
    while ~next && ishandle(hf)
        d1 = d(:,:,1) > a(1) &  d(:,:,1) < a(2);
        d2 = d(:,:,2) > b(1) &  d(:,:,2) < b(2);
        d3 = d(:,:,3) > c(1) &  d(:,:,3) < c(2);
        bw = d1 & d2 & d3;
        imshow(bw,'InitialMagnification','fit');
        drawnow limitrate        
        pause(0.1);
        ishandle(hf);
    end
    clc
    k = k + 1
end

out = [a;b;c];
end

function keyDownListener(object, eventdata)
    global a b c next shift
    key = eventdata.Key
    if strcmp(key,'q')
        a(1) = max(0,a(1)-shift*0.01);
    elseif strcmp(key,'w')
        a(1) = min(1,a(1)+shift*0.01);
    elseif strcmp(key,'a')
        b(1) = max(0,b(1)-shift*0.01);
    elseif strcmp(key,'s')
        b(1) = min(1,b(1)+shift*0.01);
    elseif strcmp(key,'z')
        c(1) = max(0,c(1)-shift*0.01);
    elseif strcmp(key,'x')
        c(1) = min(1,c(1)+shift*0.01);
    elseif strcmp(key,'i')
        a(2) = max(0,a(2)-shift*0.01);
    elseif strcmp(key,'o')
        a(2) = min(1,a(2)+shift*0.01);
    elseif strcmp(key,'k')
        b(2) = max(0,b(2)-shift*0.01);
    elseif strcmp(key,'l')
        b(2) = min(1,b(2)+shift*0.01);
    elseif strcmp(key,'comma')
        c(2) = max(0,c(2)-shift*0.01);
    elseif strcmp(key,'period')
        c(2) = min(1,c(2)+shift*0.01);
    elseif strcmp(key,'rightarrow')
        next = 1;
    elseif strcmp(key,'shift')
        shift = 5;
    end
end

function keyUpListener(object, eventdata)
    global shift
    key = eventdata.Key
    if strcmp(key,'shift')
        shift = 5;
    end
end

