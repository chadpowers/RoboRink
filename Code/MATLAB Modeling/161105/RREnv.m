classdef RREnv < handle
    %RREnv, this class simulates playing on an air hockey table,
    %can be given 3 different modes (play, watch, train) depending on
    %the use case
    
    properties (SetAccess = private)
        % PHYSICAL PROPERTIES
        puck;   % puck object
        p1;     % bottom player (human or dumb AI)
        p2;     % top player    (AI)
        table;  % tableframe
        e;      % elasticity of impact
        score;  % score
        speed;  % speed variable
        mode;   % mode
                    % 1 = play with human
                    % 2 = watch with dumb AI
                    % 3 = train
        mouse;
        bias;
        
        % DRAWING PROPERTIES
        fig;
        createMov; % whether is outputs a movie file
        
        % TIMER
        simTimer; % timer for simulation
        cycTimer; % timer for current cycle
    end
    
    properties (Access = public)
        brain;  % brain used to evaluate AI
    end
    
    methods
        
        % CONSTRUCTOR %
        function ENV = RREnv(puck,mallet,table,e,speed,mode,createMov)
            % INITIALIZE VARIABLES
            ENV.puck.d     = puck(1);
            ENV.puck.m     = puck(2);
            ENV.p1.d       = mallet(1);
            ENV.p2.d       = mallet(1);
            ENV.p1.m       = mallet(2);
            ENV.p2.m       = mallet(2);
            ENV.p1.v       = [0,0];
            ENV.p2.v       = [0,0];
            ENV.table.x    = table(1);
            ENV.table.y    = table(2);
            ENV.table.g    = table(3);
            ENV.e          = e;
            ENV.speed      = speed;
            ENV.createMov  = createMov;
            ENV.mode       = mode;
            
            ENV.score      = [0,0];
            ENV.mouse.x    = 0;
            ENV.mouse.y    = 0;
            
            
        end % END RREnv()
        
        function setUp(ENV) 
            % SETUP PLAYING FIELD
            ENV.fig = figure('units','normalized','position',[0 0 .3 .8]);
            set (ENV.fig, 'WindowButtonMotionFcn', @(object, eventdata)mouseMove(ENV,object,eventdata));
            set (ENV.fig, 'WindowButtonDownFcn', @(object, eventdata)mouseClick(ENV,object,eventdata));
            if ENV.mode == 3
                ENV.fig.Visible = 'off';
            end
            
            title(sprintf('Bottom   %d - %d   Top',ENV.score))
            axis equal
            hold on;
            ENV.table.obj = rectangle('Position',[0,0,ENV.table.x,ENV.table.y]);
            rectangle('Position',[ENV.table.x/2-ENV.table.g/2,-1,ENV.table.g,2],'FaceColor','w','EdgeColor','none');
            rectangle('Position',[ENV.table.x/2-ENV.table.g/2,ENV.table.y-1,ENV.table.g,2],'FaceColor','w','EdgeColor','none');
            ENV.puck.obj = rectangle('Curvature',[1,1],'Position',[3.9,10,ENV.puck.d,ENV.puck.d],'FaceColor','r');
            ENV.p1.obj = rectangle('Curvature',[1,1],'Position',[ENV.table.x/2,0,ENV.p1.d,ENV.p1.d],'FaceColor','b');
            ENV.p2.obj = rectangle('Curvature',[1,1],'Position',[ENV.table.x/2 - ENV.p1.d/2,ENV.table.y - ENV.p2.d - 3,ENV.p2.d,ENV.p2.d],'FaceColor','g');
            
            xlim([-1, ENV.table.x+1])
            ylim([-1, ENV.table.y+1])
        end
        
        function brain = simulate(ENV,time)
            setUp(ENV);
            
            ENV.bias = 2*(rand()-0.5);
            
            ENV.puck.obj.Position(1:2) = [ENV.table.x/2,ENV.table.y/2];
            ENV.puck.v = [0 -5];
            totalTime = 0;
            drawnow
            if ENV.mode ~= 3
                pause(1)
            end
            ENV.cycTimer = tic;
            while ishandle(ENV.fig) && (totalTime < time/ENV.speed)
                dt = ENV.speed*toc(ENV.cycTimer);
                if ENV.mode == 3
                    dt = 1/60;
                end
                ENV.cycTimer = tic;
                totalTime = totalTime + dt;
                g = step(ENV,dt);
                drawnow limitrate
                
                if ENV.mode ~= 3
                    pause(0.01); % THROTTLE BACK
                end
                
                if sum(g) == 1
                    ENV.score = ENV.score + g;
                    title(sprintf('Human   %d - %d   Computer',ENV.score))
                    pause(0.1)
                    ENV.puck.obj.Position(1:2) = [ENV.table.x/2,ENV.table.y/2];
                    ENV.puck.v = [0,-5];
                    ENV.cycTimer = tic;
                end
                
            end
            if ishandle(ENV.fig)
                close(ENV.fig)
            end
            brain = ENV.brain;
        end
        
        function [g] = step(ENV,dt)
            
            % MOVE PUCK
            ENV.puck.obj.Position(1:2) = ENV.puck.obj.Position(1:2) + dt*ENV.puck.v;
            if ENV.puck.v(1) ~= 0
                    ENV.puck.v(1) = ENV.puck.v(1)*(1-0.1*dt/abs(ENV.puck.v(1)));
            end
            if ENV.puck.v(2) ~= 0
                ENV.puck.v(2) = ENV.puck.v(2)*(1-0.1*dt/abs(ENV.puck.v(2)));
            end
            
            % MOVE P1 MALLET
            if ENV.mode == 1
                ENV.p1.newp = [max(0,min(ENV.table.x-ENV.p1.d,ENV.mouse.x-ENV.p1.d/2)),...
                    max(0,min(ENV.table.y/2-ENV.p1.d,ENV.mouse.y-ENV.p1.d/2))];
            else
                ENV.p1.newp = [max(0+1.2*ENV.puck.d,min(ENV.table.x-ENV.p1.d-1.2*ENV.puck.d,ENV.puck.obj.Position(1)+ENV.puck.d/2-ENV.p2.d/2)-ENV.bias),...
                    max(min(1*ENV.table.y/4,ENV.puck.obj.Position(2)+ENV.puck.d/2-2*rand()),0)];
            end
            ENV.p1.v = [(ENV.p1.newp(1)-ENV.p1.obj.Position(1))/dt, ...
                (ENV.p1.newp(2)-ENV.p1.obj.Position(2))/dt];
            if ENV.mode == 1
                if norm(ENV.p1.v) > 20
                    ENV.p1.v = ENV.p1.v/norm(ENV.p1.v)*20;
                end
            else
                if norm(ENV.p1.v) > 4
                    ENV.p1.v = ENV.p1.v/norm(ENV.p1.v)*4;
                end
            end
            ENV.p1.obj.Position(1:2) = ENV.p1.obj.Position(1:2)+ENV.p1.v*dt;
            
            % VALIDATE POSITION
            if ENV.p1.obj.Position(1) < 0
                ENV.p1.obj.Position(1) = 0;
                ENV.p1.v(1) = 0;
            elseif ENV.p1.obj.Position(1) > ENV.table.x - ENV.p1.d
                ENV.p1.obj.Position(1) = ENV.table.x - ENV.p1.d;
                ENV.p1.v(1) = 0;
            end
            if ENV.p1.obj.Position(2) < 0
                ENV.p1.obj.Position(2) = 0;
                ENV.p1.v(2) = 0;
            elseif ENV.p1.obj.Position(2) > ENV.table.y - ENV.p1.d
                ENV.p1.obj.Position(2) = ENV.table.y - ENV.p1.d;
                ENV.p1.v(2) = 0;
            end


            %basic AI
            input = [ENV.puck.obj.Position(1:2),ENV.puck.v,ENV.p2.obj.Position(1:2)];
            output = Genome.evaluateNet(ENV.brain.gene,input,.5);
            ENV.p2.v = 20*(output-0.5)';
            newp = ENV.p2.obj.Position(1:2)+ENV.p2.v*dt;
            if newp(1) < 0
                newp(1) = 0;
                ENV.p2.v(1) = 0;
                %ENV.brain.fitness = ENV.brain.fitness - 0.1/10000;
            elseif newp(1) > ENV.table.x-ENV.p2.d
                newp(1) = ENV.table.x-ENV.p2.d;
                ENV.p2.v(1) = 0;
                %ENV.brain.fitness = ENV.brain.fitness - 0.1/10000;
            end
            if newp(2) < ENV.table.y/2
                newp(2) = ENV.table.y/2;
                ENV.p2.v(2) = 0;
                %ENV.brain.fitness = ENV.brain.fitness - 0.1/10000;
            elseif newp(2) > ENV.table.y-ENV.p2.d
                newp(2) = ENV.table.y-ENV.p2.d;
                ENV.p2.v(2) = 0;
                %ENV.brain.fitness = ENV.brain.fitness - 0.05/10000;
            end
    
            ENV.p2.obj.Position(1:2) = newp;
    
            % detect collision
            collision(ENV);
    
            % detect goal
            g = goal(ENV);
        end
        
        function collision(ENV)
            % Collision w/ Wall
            x = ENV.puck.obj.Position(1);
            y = ENV.puck.obj.Position(2);
            if (x >= ENV.table.x - ENV.puck.d) 
                ENV.puck.v(1) = -ENV.e*ENV.puck.v(1);
                ENV.puck.obj.Position(1) = ENV.table.x - ENV.puck.d;
            elseif (x <=                0) 
                ENV.puck.v(1) = -ENV.e*ENV.puck.v(1);
                ENV.puck.obj.Position(1) = 0;
            end;
            goals = [ENV.table.x/2 - ENV.table.g/2, ENV.table.x/2 + ENV.table.g/2];
            if (y >= ENV.table.y - ENV.puck.d)
                post_dist = [x+ENV.puck.d/2-goals(1),x+ENV.puck.d/2-goals(2);y+ENV.puck.d/2-ENV.table.y,y+ENV.puck.d/2-ENV.table.y];
                [Y,I] = min([norm(post_dist(:,1)),norm(post_dist(:,2))]);
                if (x + ENV.puck.d/2 < goals(1) || x + ENV.puck.d/2 > goals(2))
                    ENV.puck.v(2) = -ENV.e*ENV.puck.v(2);
                    ENV.puck.obj.Position(2) = ENV.table.y-ENV.puck.d;
                elseif Y < ENV.puck.d/2
                    J = (1+ENV.e)*dot(ENV.puck.v,post_dist(:,I)')*post_dist(:,I)'/Y;
                    ENV.puck.v = ENV.puck.v - J;
                    m = -(ENV.puck.d/2/Y-1)*post_dist(:,I)';
                    m = [m,0,0];
                    ENV.puck.obj.Position = ENV.puck.obj.Position - m;
                end
                %ENV.brain.fitness = ENV.brain.fitness - 1/(1+pdist([ENV.puck.obj.Position(1:2)+ENV.puck.d/2;ENV.table.x/2,ENV.table.y])^2);
            end

            if (y <=                0) 
                post_dist = [x+ENV.puck.d/2-goals;y+ENV.puck.d/2,y+ENV.puck.d/2];
                [Y,I] = min([norm(post_dist(:,1)),norm(post_dist(:,2))]);
                if (x + ENV.puck.d/2 < goals(1) || x + ENV.puck.d/2 > goals(2))
                    ENV.puck.v(2) = -ENV.e*ENV.puck.v(2);
                    ENV.puck.obj.Position(2) = 0;
                elseif Y < ENV.puck.d/2
                    J = (1+ENV.e)*dot(ENV.puck.v,post_dist(:,I)')*post_dist(:,I)'/Y;
                    ENV.puck.v = ENV.puck.v - J;
                    m = -(ENV.puck.d/2/Y-1)*post_dist(:,I)';
                    m = [m,0,0];
                    ENV.puck.obj.Position = ENV.puck.obj.Position - m;
                end
                %ENV.brain.fitness = ENV.brain.fitness + 5/(1+pdist([ENV.puck.obj.Position(1:2)+ENV.puck.d/2;ENV.table.x/2,0])^2);
            end;


            % Collision w player mallet
            [ENV.puck,c] = ENV.collisionwmallet(ENV.puck,ENV.p1,ENV.e);
            if c
                %ENV.brain.fitness = ENV.brain.fitness + 0.05/pdist([ENV.puck.obj.Position(1:2)+ENV.puck.d/2;ENV.table.x/2,0]);
            end
            [ENV.puck,c] = ENV.collisionwmallet(ENV.puck,ENV.p2,ENV.e);
            if c
                %ENV.brain.fitness = ENV.brain.fitness - 2*ENV.puck.v(2)-abs(ENV.puck.v(2));
            end
        end

        


        function goal = goal(ENV)
            goal = [0,0];
            if ENV.puck.obj.Position(2) + ENV.puck.d/2 > ENV.table.y
                goal = [1,0];
                ENV.brain.fitness = ENV.brain.fitness - 15;
            elseif ENV.puck.obj.Position(2) + ENV.puck.d/2 < 0
                goal = [0,1];
                %ENV.brain.fitness = ENV.brain.fitness + 10;
            end
        end
        
        function mouseMove(ENV,object, eventdata)
            C = get (gca, 'CurrentPoint');
            ENV.mouse.x = C(1,1);
            ENV.mouse.y = C(2,2);
        end

        function mouseClick(ENV,object, eventdata)
            C = get (gca, 'CurrentPoint');
            ENV.mouse.x = C(1,1);
            ENV.mouse.y = C(2,2);

            ENV.mouse
            ENV.puck.obj.Position
            ENV.puck.v

            ENV.p1.obj.Position
            ENV.p1.v
        end
        
        function a = v(ENV) 
            a = 1;
        end
        
    end
    
    methods (Static)
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
        
    end % END METHODS
end

