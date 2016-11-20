%%%%%%
%Simulation
clear all
clc

% SIMULATION VARIABLES
puck   = [1.5, .1];
mallet = [2,5];
table  = [12 24 4];
e      = 0.9;
mode   = 3;
speed  = 1.5;
movie  = 0;



%AI VARIABLES
genomefile = 'genome_16_11_01_23.mat';
% ------ OR -------%
width    = 20;
height   = 2;
input    = 6;
output   = 2;
mutab    = 0.004;
cross    = 0.75;
pop      = 100;
elitism  = 0.05;

try
   load(genomefile)
   G
catch
    G = Genome(height,width,input,output,pop,mutab,cross,elitism)
end

if mode == 1 || mode == 2
    env = RREnv(puck,mallet,table,e,speed,mode,movie);
    env.brain = newBrain(G);
    simulate(env,Inf)
elseif mode == 3
    c = 1;
    percent = 0;
    while percent <= 75
        fprintf('\n\n')
        disp('--------------------------------------------')
        fprintf('               NEW LOOP - %d    \n',c)
        disp('--------------------------------------------')
        fprintf('\n\n')
        
        n = G.pop - size(G.newRank,1);
        tic
        
        % watch one every generation
        b = G.newBrain();
        try
            v(env(1));
            env(1).brain = b;
        catch
            env(1) = RREnv(puck,mallet,table,e,speed,2,movie);
            env(1).brain = b;
        end
        
        for i = 2:n
            b = G.newBrain();
            try
                v(env(i));
                env(i).brain = b;
            catch
                env(i) = RREnv(puck,mallet,table,e,speed,mode,movie);
                env(i).brain = b;
            end
        end
        t = sqrt(G.genCount);
        toc
        B(1) = simulate(env(1),10+t);
        toc
        disp('Visual Trial Completed...')
        disp('')
        disp('Evaluating next generation in parallel...')
        parfor i = 2:n
            B(i) = simulate(env(i),10+t);
        end
        
        for i = 1:n
            insertBrain(G,B(i))
        end
        clc
        G
        toc
        fprintf('Average time is %.2f seconds\n',toc/n)
        fprintf('The top fitness rank is %d  <->  %8f \n',G.rank(1,:))
        fprintf('The average fitness rank is %8f \n',mean(G.rank(:,2)))
        conv = convergence(G);
        fprintf('Convergence score = %g ; percent = %8.4f',conv)
        
        percent = conv(2);
        
        c = c + 1;
    end
end
    


