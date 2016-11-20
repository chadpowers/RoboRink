%%%%%%
%Follow_1D.m
clear all
clc

% SIMULATION VARIABLES
puck   = [1.5, .1];
mallet = [2,5];
table  = [12 24 4];
e      = 0.9;
mode   = 1;
speed  = 1;
movie  = 0;



%AI VARIABLES
loadfile = 'genomes/airhockey_defense1_16_11_05_15.mat';
% ------ OR -------%
width     = 6;
height    = 1;
input     = 6;
output    = 2;
mutab     = 0.001;
cross     = 0.75;
pop       = 200;
elitism   = 0.1;
writefile = 'airhockey_defense1'; %without the .mat ending, leave blank for standard conv

try
   load(loadfile)
   G
catch
    G = Genome(height,width,input,output,pop,mutab,cross,elitism,writefile)
end


if mode == 1 || mode == 2
    env = RREnv(puck,mallet,table,e,speed,mode,movie);
    env.brain = bestBrain(G);
    simulate(env,Inf)
elseif mode == 3
    c = 1;
    percent = 0;

    while percent <= 95 && G.genCount <= 100
        fprintf('\n\n')
        disp('--------------------------------------------')
        fprintf('               NEW LOOP - %d    \n',c)
        disp('--------------------------------------------')
        fprintf('\n\n')
        
        n = G.pop - size(G.newRank,1);
        tic
        
        % watch one every generation
        b = G.bestBrain();
        try
            v(env(1));
            env(1).brain = b;
        catch
            env(1) = RREnv(puck,mallet,table,e,speed,2,movie);
            env(1).brain = b;
        end
        
        t = sqrt(G.genCount);
        toc
        B(1) = simulate(env(1),10);
        toc
        
        for i = 1:n
            b = G.newBrain();
            try
                v(env(i));
                env(i).brain = b;
            catch
                env(i) = RREnv(puck,mallet,table,e,speed,mode,movie);
                env(i).brain = b;
            end
        end
        
        disp('Visual Trial Completed...')
        disp('')
        fprintf('Visual fitness score was...   %f\n',B(1).fitness)
        disp('Evaluating next generation in parallel...')
        parfor i = 1:n
            B(i) = simulate(env(i),2.5*60);
            %fprintf('%d, \n',i)
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
        hold off
        yyaxis left
        plot(G.convergenceTrend(:,1),G.convergenceTrend(:,2))
        yyaxis right
        plot(G.topRawFitnessTrend(:,1),G.topRawFitnessTrend(:,2))
        hold on
        plot(G.medRawFitnessTrend(:,1),G.medRawFitnessTrend(:,2))
        plot(G.botRawFitnessTrend(:,1),G.botRawFitnessTrend(:,2))
        
        low = max(1,size(G.topRawFitnessTrend,1)-10);
        coeffs = polyfit(G.topRawFitnessTrend(low:end,1), G.topRawFitnessTrend(low:end,2), 1);
        legend('Convergence','MAX Raw Fitness','MED Raw Fitness','MIN Raw Fitness')
        legend('Location','northwest')
        title(writefile)
        
        c = c + 1;
        if 0 == 1 %mod(c,30)
            disp('--------------------------------------------')
            fprintf('               RESTING    \n')
            disp('--------------------------------------------')
            pause(5*60)
        elseif 0 == mod(c,10)
            disp('--------------------------------------------')
            fprintf('               RESTING    \n')
            disp('--------------------------------------------')
            for i = 1:5
                fprintf('%d, ',i)
                pause(60)
            end
        end
    end
end
    


