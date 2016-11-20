function [] = play()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    draw()
    global mouse fig table puck mallet score;
    mouse.x = 0; mouse.y = 0;
    set (fig, 'WindowButtonMotionFcn', @mouseMove);
    set (fig, 'WindowButtonDownFcn', @mouseClick);
    pause(1)
    isPlaying = 1;
    tic;
    nFrames = 400;
    mov(1:nFrames) = struct('cdata', [],...
                        'colormap', []);
                    
    xlim([-1, table.x+1])
    ylim([-1, table.y+1])
    
    i = 1;
    while 1
        dt = toc;
        tic;
        g = step(dt);
        draw()
        
        pause(0.01)
        
        if sum(g) == 1
            score = score + g;
            title(sprintf('Human   %d - %d   Computer',score))
            
            pause(2)
            puck.obj.Position = [5,20,puck.d,puck.d];
            puck.v = [0,-3];
        end
        %[play,winner] = step(fig);
        
%         
%         if i < nFrames
%             mov(i) = getframe(fig);
%         elseif i == nFrames
%             mov(i) = getframe(fig);
%             vidObj = VideoWriter('simulation.mp4');
%             open(vidObj);
%             writeVideo(vidObj,mov)
%             close(vidObj);
%         end
%         i = i+1
    end
    
end

function [g] = step(dt)
    global mouse table puck mallet comp;
    % move
    puck.obj.Position = [max(0,min(puck.obj.Position(1) + dt*puck.v(1),table.x-puck.d)),...
        puck.obj.Position(2) + dt*puck.v(2),...
        puck.d,puck.d];
    if puck.v(1) ~= 0
            puck.v(1) = puck.v(1)/abs(puck.v(1))*(abs(puck.v(1))-0.1*dt);
    end
    if puck.v(2) ~= 0
        puck.v(2) = puck.v(2)/abs(puck.v(2))*(abs(puck.v(2))-0.1*dt);
    end
    
    mallet.newp = [max(0,min(table.x-mallet.d,mouse.x-mallet.d/2)),...
        max(0,min(table.y/2-mallet.d,mouse.y-mallet.d/2))];
    mallet.v = [(mallet.newp(1)-mallet.obj.Position(1))/dt, ...
        (mallet.newp(2)-mallet.obj.Position(2))/dt];
    if norm(mallet.v) > 25
        mallet.v = mallet.v/norm(mallet.v)*25;
    end
    mallet.obj.Position = [mallet.obj.Position(1:2)+mallet.v*dt, mallet.d, mallet.d];
    
    
    %basic AI
    comp.newp = [max(0,min(table.x-mallet.d,puck.obj.Position(1)+puck.d/2-comp.d/2)),...
        min(max(3*table.y/4,puck.obj.Position(2)+puck.d/2),table.y-comp.d)];
    comp.v = [(comp.newp(1)-comp.obj.Position(1))/dt, ...
        (comp.newp(2)-comp.obj.Position(2))/dt];
    if norm(comp.v) > 7
        comp.v = comp.v/norm(comp.v)*7;
    end;
    comp.obj.Position = [comp.obj.Position(1:2)+comp.v*dt, comp.d, comp.d];
    % detect collision
    collision();
    
    % detect goal
    g = goal();
end

function draw()
    drawnow limitrate
end

function collision()
    global table puck mallet goalsize comp e;
    
    % Collision w/ Wall
    x = puck.obj.Position(1);
    y = puck.obj.Position(2);
    if (x >= table.x - puck.d) puck.v(1) = -e*puck.v(1);end;
    if (x <=                0) puck.v(1) = -e*puck.v(1);end;
    goals = [table.x/2 - goalsize/2, table.x/2 + goalsize/2];
    if (y >= table.y - puck.d)
        post_dist = [x+puck.d/2-goals(1),x+puck.d/2-goals(2);y+puck.d/2-table.y,y+puck.d/2-table.y];
        [Y,I] = min([norm(post_dist(:,1)),norm(post_dist(:,2))]);
        if (x + puck.d/2 < goals(1) || x + puck.d/2 > goals(2))
            puck.v(2) = -e*puck.v(2);
            puck.obj.Position(2) = table.y-puck.d;
        elseif Y < puck.d/2
            J = (1+e)*dot(puck.v,post_dist(:,I)')*post_dist(:,I)'/Y;
            puck.v = puck.v - J;
            m = -(puck.d/2/Y-1)*post_dist(:,I)';
            m = [m,0,0];
            puck.obj.Position = puck.obj.Position - m;
        end
    end
    
    if (y <=                0) 
        post_dist = [x+puck.d/2-goals;y+puck.d/2,y+puck.d/2];
        [Y,I] = min([norm(post_dist(:,1)),norm(post_dist(:,2))]);
        if (x + puck.d/2 < goals(1) || x + puck.d/2 > goals(2))
            puck.v(2) = -e*puck.v(2);
            puck.obj.Position(2) = 0;
        elseif Y < puck.d/2
            J = (1+e)*dot(puck.v,post_dist(:,I)')*post_dist(:,I)'/Y;
            puck.v = puck.v - J;
            m = -(puck.d/2/Y-1)*post_dist(:,I)';
            m = [m,0,0];
            puck.obj.Position = puck.obj.Position - m;
        end
    end;
    
    
    % Collision w player mallet
    puck = collisionwmallet(puck,mallet,e);
    if 0 %(pdist([puck.obj.Position(1)+puck.d/2,puck.obj.Position(2)+puck.d/2;...
            %mallet.obj.Position(1)+mallet.d/2,mallet.obj.Position(2)+mallet.d/2])<=(mallet.d+puck.d)/2)
        point = ((mallet.obj.Position(1:2)+mallet.d/2)*puck.d + (puck.obj.Position(1:2)+puck.d/2)*mallet.d)/(mallet.d+puck.d);
        
        mov = point - [puck.obj.Position(1)+puck.d/2,puck.obj.Position(2)+puck.d/2];
        mov_n = mov/norm(mov);
        J = (1+e)*mallet.m/(mallet.m+puck.m)*dot(puck.v-mallet.v,mov_n)*mov_n;
        
        m = -(puck.d/2/norm(mov)-1)*mov;
        m = [m,0,0];
        puck.obj.Position = puck.obj.Position + m;
        
        puck.v = puck.v - J;
        %puck.v = (puck.v*(1-10)+2*10*mallet.v)/(1+10);
        
        %puck.v = 2*mallet.v-puck.v;
        
    end
    puck = collisionwmallet(puck,comp,e);
end

function [p] = collisionwmallet(puck, mallet, e)
    if (pdist([puck.obj.Position(1)+puck.d/2,puck.obj.Position(2)+puck.d/2;...
            mallet.obj.Position(1)+mallet.d/2,mallet.obj.Position(2)+mallet.d/2])<=(mallet.d+puck.d)/2)
        point = ((mallet.obj.Position(1:2)+mallet.d/2)*puck.d + (puck.obj.Position(1:2)+puck.d/2)*mallet.d)/(mallet.d+puck.d);
        
        mov = point - [puck.obj.Position(1)+puck.d/2,puck.obj.Position(2)+puck.d/2];
        mov_n = mov/norm(mov);

        J = (1+e)*mallet.m/(mallet.m+puck.m)*dot(puck.v-mallet.v,mov_n)*mov_n;
        
        m = -(puck.d/2/norm(mov)-1)*mov;
        m = [m,0,0];
        puck.obj.Position = puck.obj.Position + m;
        
        puck.v = puck.v - J;
    end
    p = puck;
end


function [goal] = goal()
    global puck table;
    goal = [0,0];
    if puck.obj.Position(2) + puck.d/2 > table.y
        goal = [1,0];
    elseif puck.obj.Position(2) + puck.d/2 < 0
        goal = [0,1];
    end
end

function mouseMove(object, eventdata)
    global mouse;
    C = get (gca, 'CurrentPoint');
    mouse.x = C(1,1);
    mouse.y = C(2,2);
end

function mouseClick(object, eventdata)
    global mouse puck mallet table fig;
    C = get (gca, 'CurrentPoint');
    mouse.x = C(1,1);
    mouse.y = C(2,2);
    
    mouse
    puck.obj.Position
    puck.v
    
    mallet.obj.Position
    mallet.v
end