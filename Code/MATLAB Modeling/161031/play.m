function [] = play()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    draw()
    global mouse fig table puck score mallet;
    global width height input output mutability crossover population;
    global comp genome rank;
    
    mouse.x = 0; mouse.y = 0;
    set (fig, 'WindowButtonMotionFcn', @mouseMove);
    set (fig, 'WindowButtonDownFcn', @mouseClick);
    pause(1)
    isPlaying = 1;
    tic;
    
    compTime = tic;
    genome = createGenome(width,height,input,output,population);
    rank = [[1:population]',zeros(population,1)-100];
    comp.brain = genome(1);
    
    nFrames = 400;
    mov(1:nFrames) = struct('cdata', [],...
                        'colormap', []);
                    
    xlim([-1, table.x+1])
    ylim([-1, table.y+1])
    
    i = 1;
    while 1
        dt = toc;
        tic;
        if toc(compTime) > 10
            %update AI
            comp.brain = newBrain();
            puck.obj.Position = [5,14,puck.d,puck.d];
            puck.v = [0,-7];
            mallet.obj.Position(1:2) = [table.x/2,0];
            comp.obj.Position(1:2) = [table.x/2,table.y];
            rank;
            save genome.mat genome
            % reset Timer
            compTime = tic;
            tic;
        end
        
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
    
    
    %mallet.newp = [max(0,min(table.x-mallet.d,mouse.x-mallet.d/2)),...
    %    max(0,min(table.y/2-mallet.d,mouse.y-mallet.d/2))];
    mallet.newp = [max(0,min(table.x-mallet.d,puck.obj.Position(1)+puck.d/2-comp.d/2)-2*rand()),...
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
    comp.v = 7*evaluateNet(comp.brain.gene,input)';
    %comp.v = [(comp.newp(1)-comp.obj.Position(1))/dt, ...
        %(comp.newp(2)-comp.obj.Position(2))/dt];
    if norm(comp.v) > 7
        comp.v = comp.v/norm(comp.v)*7;
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
        comp.brain.fitness = comp.brain.fitness - 50;
    elseif puck.obj.Position(2) + puck.d/2 < 0
        goal = [0,1];
        comp.brain.fitness = comp.brain.fitness + 50;
    end
end

function [genome] = createGenome(width,height,input,output,population)
    genome(population).gene(height+1) = 0;
    for i = 1:population
        gene(height+1).m = 0;
        if height > 0
            gene(1).m = 2*rand(width,input+1)-1;
            j = 2;
            while j <= height
                gene(j).m = 2*rand(width,width+1)-1;
                j = j+1;
            end
            gene(j).m = 2*rand(output,width+1)-1;
        else
            gene(1).m = 2*rand(output,input+1)-1;
        end
        genome(i).gene = gene;
        genome(i).fitness = 0;
        genome(i).place = i;
    end
end

function [output] = evaluateNet(net,input)
    output = input';
    for i = 1:size(net,2)
        output = net(i).m*[output;-1];
    end       
end

function [brain] = newBrain()
    global height width input output;
    global genome rank comp mutability crossover population;
    
    % insert latest genome
    if comp.brain.fitness > rank(end,2)
        i = 1;
        while comp.brain.fitness < rank(i,2)
            i = i+1;
        end
        if i == 1
            rank = [comp.brain.place, comp.brain.fitness; rank(i:end-1,:)];
        elseif i < population + 1
            rank = [rank(1:i-1,:); comp.brain.place, comp.brain.fitness; rank(i:end-1,:)];
        end
    end
    
    worst = rank(end,1);
    
    % pick 2 new genomes
    n = size(rank,1)*(size(rank,1)-1)/2*rand();
    
    
    p1 = 1;
    t = population;
    while t < n
        p1 = p1 + 1;
        t = t + 40 - p1;
    end
    
    n = size(rank,1)*(size(rank,1)-1)/2*rand();
    p2 = 1;
    t = population;
    while t < n
        p2 = p2 + 1;
        t = t + 40 - p2;
    end
    
    % cross
    gene(height+1).m = 0;
    if height > 0
        rows = rand(width,1);
        gene(1).m = zeros(width,input+1);
        gene(1).m(rows<crossover,:)=genome(p1).gene(1).m(rows<crossover,:);
        gene(1).m(rows>=crossover,:)=genome(p2).gene(1).m(rows>=crossover,:);
        j = 2;
        while j <= height
            rows = rand(width,1);
            gene(j).m = zeros(width,width+1);
            gene(j).m(rows<crossover,:)=genome(p1).gene(j).m(rows<crossover,:);
            gene(j).m(rows>=crossover,:)=genome(p2).gene(j).m(rows>=crossover,:);
            j = j+1;
        end
        rows = rand(output,1);
        gene(j).m = zeros(output,width+1);
        gene(j).m(rows<crossover,:)=genome(p1).gene(j).m(rows<crossover,:);
        gene(j).m(rows>=crossover,:)=genome(p2).gene(j).m(rows>=crossover,:);
    else
        rows = rand(output,1);
        gene(1).m = zeros(output,input+1);
        gene(1).m(rows<crossover,:)=genome(p1).gene(1).m(rows<crossover,:);
        gene(1).m(rows>=crossover,:)=genome(p2).gene(1).m(rows>=crossover,:);
    end
    
    % mutate
    for i = 1:height+1
        for j = size(gene(i).m,1)
            for k = size(gene(i).m,2)
                if rand() < mutability
                    gene(i).m(j,k) = 2*rand()-1;
                end
            end
        end
    end
            
    genome(worst).gene = gene;
    genome(worst).fitness = 0;
    genome(worst).place = worst;
    
    brain = genome(worst);
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
