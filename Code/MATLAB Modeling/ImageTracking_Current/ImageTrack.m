% First pass at image tacking
% ImageTrack.m

swingObj = VideoReader('Vd1.mov');


%Then, determine the width and height of the video.
vidWidth = swingObj.Width;
vidHeight = swingObj.Height;
scale = 1/10;

try 
    load('rawmov.mat')
catch

    %Create a movie structure array, mov.
    mov = struct('cdata',zeros(vidHeight*scale,vidWidth*scale,3,'uint8'),...
        'colormap',[]);
    
    
    k = 1;
    while hasFrame(swingObj)
        values = readFrame(swingObj);
        d = 255*rgb2hsv(values);
        d = imresize(d,scale);
        mov(k).cdata = uint8(d);
        imshow(mov(k).cdata);
        drawnow limitrate
        clc
        k = k + 1
    end
    
    save('rawmov.mat','mov')
end

bwmov = struct('cdata',zeros(vidHeight*scale,vidWidth*scale,1,'uint8'),...
        'colormap',[]);

%%
centroids = zeros(size(mov,2),2);
%Read one frame at a time until the end of the video is reached.
tic
for k = 1:size(mov,2)
    d = int16(mov(k).cdata);
    clc
    k
    
    % move red away from boundaryc
    d(:,:,1) = d(:,:,1) - 50;
    d(d(:,:,1)<0) = d(d(:,:,1)<0) + 256;
    dis1 = d(:,:,1) - 201;
    dis2 = d(:,:,2) - 227;
    dis3 = d(:,:,3) - 202;
    
    
    v = -double((max(dis1.^2+dis2/20.^2+dis3/20.^2 - 200,0))).^(1/10);
    %pause(0.05)
    %v = v.^0.1;
    
    v = v-min(min(v));
    v(:,:) = v/max(max(v));
    bw_gpu = v>0.95;
    bw = bw_gpu;%gather(bw_gpu);
    s = regionprops(bw,'centroid','area');
    if size(s,1) > 0
        [~,idx]=max([s.Area]);
        s=s(idx);
        centroids(k,:) = [s(1).Centroid];
    else
        centroids(k,:) = [0,0];
    end
%      hf = figure(3);
%      set(hf,'position',[10 10 vidWidth/2 vidHeight/2]);
%      imshow(bw)
%     hold all
%     plot(centroids(:,1), centroids(:,2),'bo');
%     plot(centroids(end,1),centroids(end,2),'gO');
%     hold off
%     drawnow
%     bwmov(k).cdata(:,:,1) = uint8((v>0.75)*255);
%     bwmov(k).cdata(:,:,2) = bwmov(k).cdata(:,:,1);
%     bwmov(k).cdata(:,:,3) = bwmov(k).cdata(:,:,1);
end

ratio = 60*1.8/(max(centroids(:,1))-min(centroids(:,1)));
a = diff(centroids)*ratio;
b = sqrt(a(:,1).^2+a(:,2).^2);

c = tsmovavg(b','s',3);

toc

plot(b)
%mov(k).cdata = mov(k).cdata


% %%
% % %Size a figure based on the width and height of the video. 
% % %Then, play back the movie once at the video frame rate.
% % hf = figure(1);
% % hf2 = figure(2);
% % 
%  startFrame = max(0,200);
%  endFrame   = min(size(mov,2),500);
%  rate       = 5;
% % 
% % playMov = mov(startFrame:endFrame);
% % 
% % set(hf,'position',[10 10 vidWidth/2 vidHeight/2]);
% % set(hf2,'position',[10 10 200 200]);
% % movie(hf,bwmov,1,swingObj.FrameRate/1);
% % 
% % % prompt1 = 'Do you want to start measuring bat/ball locations? Y/N [Y]: ';
% % %         str = input(prompt1,'s');
% % %         if isempty(str)
% % %         str = 'Y';
% % %         end
% %         
% i = 0;
% j = 1;
% hsv_values = zeros(floor((startFrame-endFrame)/rate),3);
% while i < endFrame %strcmp(str,'Y')==1
%     hf;
%     hold on
%     figure(1)
%     %set(hf,'position',[10 10 vidWidth/2 vidHeight/2]);
%     imagesc(mov(i+startFrame).cdata)
%     %axis([xi xa yi ya]);
%     ax=gca;
%     ax.YDir='rev';
%     %title(sprintf('Chad 10.14.2016, clip %d, frame %d: Press "e" for end of bat, "k" for knob of bat, "b" for ball, and "x" for no ball',c,i+start))
%     [y,x,key] = ginput(1); %the 'ob' here tells ginput you are going to track ob objects
%     x = round(x)
%     y = round(y)
%     hsv_values(j,:) = mov(i+startFrame).cdata(x,y,:);
%     hsv_values(j,:)
%     
%     figure(2)
%     clf
%     ax=gca;
%     ax.YDir='rev';
%     %imagesc(mov(i+startFrame).cdata(max(1,x-100):min(x+100,vidHeight),max(1,y-100):min(vidWidth,y+100),:))
%     %rectangle('Position',[0,0,100,100],'FaceColor',hsv_values(j,:)/256,'EdgeColor','none');
%     
%     %hold off
%     %D((ob*i-2):ob*i,:) = cat(2,x,y,key);
%     i = i + rate;
%     j = j + 1;
% %         %prompt2 = 'Do you want to measure another frame? Y/N [Y]: ';
% %         str = input(prompt2,'s');
% %         if isempty(str)
% %         str = 'Y';
% %         end
% end

