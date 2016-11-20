function [] = play()

    % set up
    global mouse fig table puck score mallet;
    global width height input output mutability crossover population;
    global comp rank;
    global bias;
    
    elitism = 0.0;
    
    bias = 0;
    
    mouse.x = 0; mouse.y = 0;
    set (fig, 'WindowButtonMotionFcn', @mouseMove);
    set (fig, 'WindowButtonDownFcn', @mouseClick);
    pause(1)
    isPlaying = 1;
    tic;
    
    compTime = tic;
    genome = Genome(height,width,input,output,population,mutability,crossover,elitism);
    rank = [[1:population]',zeros(population,1)-100];
    comp.brain = newBrain(genome);
                    
    xlim([-1, table.x+1])
    ylim([-1, table.y+1])
    
    i = 1;
    while 1
        dt = 2*toc;
        tic;
        if toc(compTime) > 3 + genome.genCount
            %update AI
            clc
            insertBrain(genome,comp.brain);
            comp.brain = newBrain(genome);
            puck.obj.Position = [5,14,puck.d,puck.d];
            puck.v = [0,-7];
            
            %genome.rank(1:end,:)
            genome.newRank(1:end,:)
            mallet.obj.Position(1:2) = [table.x/2,0];
            comp.obj.Position(1:2) = [table.x/2,table.y];
            
            bias = 3*(rand()-0.5);
            compTime = tic;
            tic;
        end
        
        g = step(dt);
        draw()
        
        pause(0.01)
        
        if sum(g) == 1
            score = score + g;
            title(sprintf('Human   %d - %d   Computer',score))
            
            pause(0.1)
            puck.obj.Position = [5,20,puck.d,puck.d];
            puck.v = [0,-3];
        end
        
    end
    
end

function [g] = step(dt)
    global mouse table puck mallet comp;
    global bias;
    global genome;
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
    
    
    %mallet.newp = [max(0,min(table.x-mallet.d,mouse.x-mallet.d/2)),...
    %    max(0,min(table.y/2-mallet.d,mouse.y-mallet.d/2))];
    mallet.newp = [max(0+1.2*puck.d,min(table.x-mallet.d-1.2*puck.d,puck.obj.Position(1)+puck.d/2-comp.d/2)-bias),...
        max(min(1*table.y/4,puck.obj.Position(2)+puck.d/2-4*rand()),0)];
    mallet.oldv = mallet.v;
    mallet.v = [(mallet.newp(1)-mallet.obj.Position(1))/dt, ...
        (mallet.newp(2)-mallet.obj.Position(2))/dt];
    %mallet.a = (mallet.v - mallet.oldv)/dt;
    if norm(mallet.v) > 4
        mallet.v = mallet.v/norm(mallet.v)*4;
    end
    %mallet.v = mallet.oldv + mallet.a*dt;
    mallet.obj.Position = [mallet.obj.Position(1:2)+mallet.v*dt, mallet.d, mallet.d];
    if mallet.obj.Position(1) < 0
        mallet.obj.Position(1) = 0;
        mallet.v(1) = 0;
    elseif mallet.obj.Position(1) > table.x - mallet.d
        mallet.obj.Position(1) = table.x - mallet.d;
        mallet.v(1) = 0;
    end
    if mallet.obj.Position(2) < 0
        mallet.obj.Position(2) = 0;
        mallet.v(2) = 0;
    elseif mallet.obj.Position(2) > table.y - mallet.d
        mallet.obj.Position(2) = table.y - mallet.d;
        mallet.v(2) = 0;
    end
    
    
    %basic AI
    comp.newp = [max(0,min(table.x-mallet.d,puck.obj.Position(1)+puck.d/2-comp.d/2)),...
        min(max(3*table.y/4,puck.obj.Position(2)+puck.d/2),table.y-comp.d)];
    input = [puck.obj.Position(1:2),puck.v,comp.obj.Position(1:2)];
    comp.v = 20*(Genome.evaluateNet(comp.brain.gene,input,.5)-0.5)';
    %comp.v = [(comp.newp(1)-comp.obj.Position(1))/dt, ...
        %(comp.newp(2)-comp.obj.Position(2))/dt];
    if norm(comp.v) > 10
        comp.v = comp.v/norm(comp.v)*10;
    end;
    newp = comp.obj.Position(1:2)+comp.v*dt;
    if newp(1) < 0;
        newp(1) = 0;
        comp.v(1) = 0;
        comp.brain.fitness = comp.brain.fitness - 0.1;
    elseif newp(1) > table.x-comp.d
        newp(1) = table.x-comp.d;
        comp.v(1) = 0;
        comp.brain.fitness = comp.brain.fitness - 0.1;
    end
    if newp(2) < table.y/2;
        newp(2) = table.y/2;
        comp.v(2) = 0;
        comp.brain.fitness = comp.brain.fitness - 0.1;
    elseif newp(2) > table.y-comp.d
        newp(2) = table.y-comp.d;
        comp.v(2) = 0;
        comp.brain.fitness = comp.brain.fitness - 0.1;
    end
    
    comp.obj.Position = [newp, comp.d, comp.d];
    
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
        comp.brain.fitness = comp.brain.fitness - 1/pdist([puck.obj.Position(1:2)+puck.d/2;table.x/2,table.y]);
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
        comp.brain.fitness = comp.brain.fitness + 1/pdist([puck.obj.Position(1:2)+puck.d/2;table.x/2,0]);
    end;
    
    
    % Collision w player mallet
    [puck,c] = collisionwmallet(puck,mallet,e);
    if c
        comp.brain.fitness = comp.brain.fitness + 1/pdist([puck.obj.Position(1:2)+puck.d/2;table.x/2,0]);
    end
    [puck,c] = collisionwmallet(puck,comp,e);
    if c
        comp.brain.fitness = comp.brain.fitness - puck.v(2);
    end
end

function [p,c] = collisionwmallet(puck, mallet, e)
    c = 0;
    potential = (puck.obj.Position(2) +   puck.d > mallet.obj.Position(2)) && ...
                (puck.obj.Position(2) - mallet.d < mallet.obj.Position(2)) && ...
                (puck.obj.Position(1) +   puck.d > mallet.obj.Position(1)) && ...
                (puck.obj.Position(1) - mallet.d < mallet.obj.Position(1));
    if (potential)
        if (pdist([puck.obj.Position(1)+puck.d/2,puck.obj.Position(2)+puck.d/2;...
            mallet.obj.Position(1)+mallet.d/2,mallet.obj.Position(2)+mallet.d/2])<=(mallet.d+puck.d)/2)
            
            c = 1;
            point = ((mallet.obj.Position(1:2)+mallet.d/2)*puck.d + (puck.obj.Position(1:2)+puck.d/2)*mallet.d)/(mallet.d+puck.d);
            
            mov = point - [puck.obj.Position(1)+puck.d/2,puck.obj.Position(2)+puck.d/2];
            mov_n = mov/norm(mov);
            
            J = (1+e)*mallet.m/(mallet.m+puck.m)*dot(puck.v-mallet.v,mov_n)*mov_n;
        
            m = -(puck.d/2/norm(mov)-1)*mov;
            m = [m,0,0];
            puck.obj.Position = puck.obj.Position + m;
        
            puck.v = puck.v - J;
        end
    end
    p = puck;
end


function [goal] = goal()
    global puck table comp;
    goal = [0,0];
    if puck.obj.Position(2) + puck.d/2 > table.y
        goal = [1,0];
        comp.brain.fitness = comp.brain.fitness - 30;
    elseif puck.obj.Position(2) + puck.d/2 < 0
        goal = [0,1];
        comp.brain.fitness = comp.brain.fitness + 10;
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
