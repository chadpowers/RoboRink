%%%%%%
%Follow_1D.m
clear all
clc

% SIMULATION VARIABLES
puck   = [1.5, .1];
mallet = [2,5];
table  = [24 6 4];
e      = 0.9;
mode   = 3;
speed  = 1.5;
movie  = 0;



%AI VARIABLES
loadfile = 'genomes/test12.mat';
% ------ OR -------%
height    = 1;
width     = 6;
input     = 6;
output    = 1;
mutab     = 0.01;
cross     = 0.6;
pop       = 200;
elitism   = 0.2;
writefile = 'test13'; %without the .mat ending, leave blank for standard conv

try
   load(loadfile)
   G
catch
    G = Genome(height,width,input,output,pop,mutab,cross,elitism,writefile)
end


if mode == 1 || mode == 2
    env = SimpleEnv(puck,mallet,table,e,speed,mode,movie);
    env.brain = newBrain(G);
    simulate(env,Inf)
elseif mode == 3
    c = 1;
    percent = 0;

    while percent <= 99
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
            env(1) = SimpleEnv(puck,mallet,table,e,speed,2,movie);
            env(1).brain = b;
        end
        
        for i = 2:n
            b = G.newBrain();
            try
                v(env(i));
                env(i).brain = b;
            catch
                env(i) = SimpleEnv(puck,mallet,table,e,speed,mode,movie);
                env(i).brain = b;
            end
        end
        t = sqrt(G.genCount);
        toc
        B(1) = simulate(env(1),20);
        toc
        disp('Visual Trial Completed...')
        disp('')
        fprintf('Visual fitness score was...   %f\n',B(1).fitness)
        disp('Evaluating next generation in parallel...')
        parfor i = 2:n
            B(i) = simulate(env(i),25);
        end
        
        for i = 1:n
            insertBrain(G,B(i))
        end
        clc
        G
        toc
        fprintf('Average time is %.2f seconds\n',toc/n)
        fprintf('The top fitness score is %d  <->  %8f \n',G.rank(1,:))
        fprintf('The average of fitness scores is %8f \n',mean(G.rank(:,2)))
        fprintf('The stddev. of fitness scores is %8f \n',std(G.rank(:,2)))
        conv = convergence(G);
        fprintf('Convergence score = %g ; percent = %8.4f',conv)
        
        percent = conv(2);
        
        figure(2)
        yyaxis left
        plot(G.convergenceTrend(:,1),G.convergenceTrend(:,2))
        yyaxis right
        plot(G.topRawFitnessTrend(:,1),G.topRawFitnessTrend(:,2))
        low = max(1,size(G.topRawFitnessTrend,1)-10);
        coeffs = polyfit(G.topRawFitnessTrend(low:end,1), G.topRawFitnessTrend(low:end,2), 1)
        legend('Convergence','Raw Fitness Score')
        legend('Location','northwest')
        title(writefile)
        
        c = c + 1;
    end
end
    


