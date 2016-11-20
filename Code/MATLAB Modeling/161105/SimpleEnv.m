classdef SimpleEnv < handle
    %RREnv, this class simulates playing on an air hockey table,
    %can be given 3 different modes (play, watch, train) depending on
    %the use case
    
    properties (SetAccess = private)
        % PHYSICAL PROPERTIES
        puck;   % puck object
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
        function ENV = SimpleEnv(puck,mallet,table,e,speed,mode,createMov)
            % INITIALIZE VARIABLES
            ENV.puck.d     = puck(1);
            ENV.puck.m     = puck(2);
            ENV.p2.d       = mallet(1);
            ENV.p2.m       = mallet(2);
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
            if ENV.mode == 3
                ENV.fig.Visible = 'off';
            end
            
            set (ENV.fig, 'WindowButtonMotionFcn', @(object, eventdata)mouseMove(ENV,object,eventdata));
            set (ENV.fig, 'WindowButtonDownFcn', @(object, eventdata)mouseClick(ENV,object,eventdata));
            
            title(sprintf('Bottom   %d - %d   Top',ENV.score))
            axis equal
            hold on;
            ENV.table.obj = rectangle('Position',[0,0,ENV.table.x,ENV.table.y]);
            ENV.puck.obj = rectangle('Curvature',[1,1],'Position',[0,ENV.table.y/2,ENV.puck.d,ENV.puck.d],'FaceColor','r');
            ENV.p2.obj = rectangle('Curvature',[1,1],'Position',[ENV.table.x/2 - ENV.p2.d/2,ENV.table.y - ENV.p2.d - 2,ENV.p2.d,ENV.p2.d],'FaceColor','g');
            
            xlim([-1, ENV.table.x+1])
            ylim([-1, ENV.table.y+1])
        end
        
        function brain = simulate(ENV,time)
            setUp(ENV);
            
            ENV.bias = 2*(rand()-0.5);
            
            ENV.puck.obj.Position(1:2) = [3,ENV.table.y/2];
            ENV.puck.v = [0 0];
            totalTime = 0;
            drawnow
            if ENV.mode ~= 3
                pause(1)
            end
            ENV.cycTimer = tic;
            while ishandle(ENV.fig) && (totalTime < time)
                dt = ENV.speed*toc(ENV.cycTimer);
                if ENV.mode == 3
                    dt = 1/60;
                end
                ENV.cycTimer = tic;
                totalTime = totalTime + dt;
                ENV.puck.v = [5*sin(totalTime/2 + ENV.bias), 0];
                step(ENV,dt);
                drawnow limitrate
                
                if ENV.mode ~= 3
                    pause(0.01); % THROTTLE BACK
                end
                
            end
            if ishandle(ENV.fig)
                close(ENV.fig)
            end
            brain = ENV.brain;
        end
        
        function step(ENV,dt)
            
            % MOVE PUCK
            ENV.puck.obj.Position(1:2) = ENV.puck.obj.Position(1:2) + dt*ENV.puck.v;


            %basic AI
            input = [ENV.puck.obj.Position(1:2),ENV.puck.v,ENV.p2.obj.Position(1:2)];
            output = (Genome.evaluateNet(ENV.brain.gene,input,0.8))';
            if size(output,2) == 1
                ENV.p2.v = [40*(output-0.5),0];
            else
                ENV.p2.v = [20*(output(1)-output(2)),0];
            end
            newp = ENV.p2.obj.Position(1:2)+ENV.p2.v*dt;
            
            ENV.brain.fitness = ENV.brain.fitness + dt/100/(1+abs(ENV.puck.obj.Position(1) - newp(1)));
    
            ENV.p2.obj.Position(1:2) = newp;
    
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

            %ENV.mouse
            ENV.puck.obj.Position
            ENV.puck.v

            ENV.p2.obj.Position
            ENV.p2.v
        end
        
        function a = v(ENV) 
            a = 1;
        end
        
    end
    
    methods (Static)
        
        
    end % END METHODS
end

