% First pass at image tacking
% ImageTrack.m
clear all 

swingObj = VideoReader('Vd2.mov');


%Then, determine the width and height of the video.
vidWidth = swingObj.Width;
vidHeight = swingObj.Height;
scale = 1/10;

try 
    load('rawmov2.mat')
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
    
    save('rawmov2.mat','mov')
end

bwmov = struct('cdata',zeros(vidHeight*scale,vidWidth*scale,1,'uint8'),...
        'colormap',[]);
    
xmin = 1000;
xmax = 0;
ymin = 1000;
ymax = 0;

%%
clear centroids d dis1 dis2 dis3 v bw bw_gpu
centroids = zeros(size(mov,2),2);


hf = figure(3);
set(hf,'position',[10 10 vidWidth/2 vidHeight/2]);
%Read one frame at a time until the end of the video is reached.
tic
for k = 500:size(mov,2)
    % PRINT LOOP NUMBER
    clc
    k
    
    m_pos = [50,53];
    
    % READ IN RGB VIDEO FRAME   AND
    % PROCESS TO BINARY
    mode = 2; % DEFINITELY USE 2, 27.7% Faster than 1
    if mode == 1
        d = int16(mov(k).cdata);
        d(:,:,1) = d(:,:,1) - 50;
        d(d(:,:,1)<0) = d(d(:,:,1)<0) + 256;
        dis1 = d(:,:,1) - 201;
        dis2 = d(:,:,2) - 227;
        dis3 = d(:,:,3) - 202;
        v = -double((max(dis1.^2+dis2/20.^2+dis3/20.^2 - 200,0))).^(1/10);
        v = v-min(min(v));
        v(:,:) = v/max(max(v));
        bw_gpu = v>0.95;
        bw = bw_gpu;
    elseif mode == 2
        d = double(mov(k).cdata)/255;
        a = [0.58, 0.93];
        b = [0.34, 1.00];
        c = [0.29, 1.00];
        d(:,:,1) = d(:,:,1) - 0.1;
        d(d(:,:,1)<0) = d(d(:,:,1)<0) + 1;
        d1 = d(:,:,1) > a(1) &  d(:,:,1) < a(2);
        d2 = d(:,:,2) > b(1) &  d(:,:,2) < b(2);
        d3 = d(:,:,3) > c(1) &  d(:,:,3) < c(2);
        bw = d1 & d2 & d3;
    end
    
    % FIND BLOB
    s = regionprops(bw,'centroid','area');
    if size(s,1) > 0
        [~,idx]=max([s.Area]);
        s=s(idx);
        centroids(k,:) = [s(1).Centroid];
    else
        centroids(k,:) = [0,0];
    end
    
    % UPDATE SIDES
    if (centroids(k,1) > xmax); xmax = centroids(k,1);
    elseif (centroids(k,1) < xmin); xmin = centroids(k,1); end;
    if (centroids(k,2) > ymax); ymax = centroids(k,2);
    elseif (centroids(k,2) < ymin); ymin = centroids(k,2); end;
    
    % VELOCITY + TRAJECTORY PREDICTION
    length = k;%size(centroids,1);
    if length > 1
        velo = mean(diff(centroids(max(1,length-3):length,:)));
        if velo(1) > 0
            time = 40;%max(0,min(5,(xmax-centroids(k,1))/velo(1)));
        else
            time = max(0,min(50,(xmin-centroids(k,1))/velo(1)));
        end
        time_variables = 0:1/60:round(time*1.1);
        trajectory = centroids(k,:)'+ velo'*time_variables;
        
        % Y BOUNCE
        l = 0;
        while l < 5 && k > 50 && size(trajectory,2) > 0 && (max(trajectory(2,:))>ymax || (min(trajectory(2,:))<ymin))
            trajectory(2,trajectory(2,:) > ymax) = 2*ymax - trajectory(2,trajectory(2,:) > ymax);
            trajectory(2,trajectory(2,:) < ymin) = 2*ymin - trajectory(2,trajectory(2,:) < ymin);
            l = l + 1;
        end
        
        % X BOUNCE
        l = 5;
        while l < 5 && k > 50 && size(trajectory,2) > 0 && (max(trajectory(1,:))>xmax || (min(trajectory(1,:))<xmin))
            trajectory(1,trajectory(1,:) > xmax) = 2*xmax - trajectory(1,trajectory(1,:) > xmax);
            trajectory(1,trajectory(1,:) < xmin) = 2*xmin - trajectory(1,trajectory(1,:) < xmin);
            l = l + 1;
        end
        
        % LOCATION DECISON
        
        if velo(1) < 0
            % binary search for lowest x value that is greater than say 50;
            % use that index's y value to move to that point
            j = 1;
            while j < size(trajectory,2) && trajectory(1,j) > 50
                j = j + 1;
            end
            
            if j == size(trajectory,2)
                some_y = 53;
            else
                some_y = trajectory(2,j);
            end
            m_pos = [50,some_y];
        else
            % move around the middle
            m_pos = [50,53];
        end
    end
    
    
    
    
    % DRAW
    draw = 1;
    if draw
        %clf
        %title(k/60)
        imshow(bw,'InitialMagnification','fit')
        hold all
        G = linspace(10,100,30)';
        myGmap = horzcat(zeros(size(G)),G/100, zeros(size(G)));
        scatter(centroids(k-29:k,1), centroids(k-29:k,2),G+.1,myGmap,'o','filled');
        scatter(m_pos(1),m_pos(2),1000,[.3,.6,1],'o','filled');
        plot(centroids(end,1),centroids(end,2),'gO');
        plot(trajectory(1,:),trajectory(2,:),'r')
        title(round(k/.6)/100)
        hold off
        drawnow
        %pause(0.1)
    end
end

tablelength = 144.6; %(max(centroids(:,1))-min(centroids(:,1)))
ratio = 60*1.8/tablelength;

clear a b c 
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

