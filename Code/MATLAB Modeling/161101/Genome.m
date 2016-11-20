classdef Genome < handle
    %GENOME, this class contains the entire gene pool for the neural net
    % AI that controls the computer mallet, genome can be shared across
    % various, simultaneous simulations of mallet
    
    properties (SetAccess = private)
        % NETWORK PROPERTIES
        height;  % height of hidden layers
        width;   % width of hidden layers
        input;   % number of inputs
        output;  % number of outputs
        pop;     % population of genome
        
        % EVOLUTION PROPERTIES
        rank;       % directory of generation based on fitness
        mutability; % percent chance of complete mutation - 10x for minor mutation
        crossover;  % percent of 1st parent that crosses over
        elitism;    % percent of previous generation that directly goes
        newRank; % directory of new generation based on fitness
        
        % HELPER/COUNTERS
        genCount; % number of generations
        brainCount;      % number of brains
    end
    
    properties (Access = private)
        nextBrain;
        curGen;  % all genes of current generation
        newGen;  % all genes of old generation
    end
    
    methods
        
        % CONSTRUCTOR %
        function G = Genome(h,w,in,out,p,mut,cross,elit)
            % INITIALIZE VARIABLES
            G.height     = h;
            G.width      = w;
            G.input      = in;
            G.output     = out;
            G.pop        = p;
            G.mutability = mut;
            G.crossover  = cross;
            G.elitism    = elit;
            G.genCount   = 1;
            G.brainCount = 0;
            G.nextBrain  = 1;
            
            % SET UP RANK
            G.rank = [[1:G.pop]',zeros(G.pop,1)-1000];
            
            % INITIALIZE POPULATION
            G.curGen(G.pop).gene(G.height+1) = 0;
            for i = 1:G.pop
                gene(G.height+1).m = 0;
                if G.height > 0
                    gene(1).m = 2*rand(G.width,G.input+1)-1;
                    j = 2;
                    while j <= G.height
                        gene(j).m = 2*rand(G.width,G.width+1)-1;
                        j = j+1;
                    end
                    gene(j).m = 2*rand(G.output,G.width+1)-1;
                else
                    gene(1).m = 2*rand(G.output,G.input+1)-1;
                end
                G.curGen(i).gene = gene;
                G.curGen(i).fitness = 0;
                G.curGen(i).place = i;
            end
        end % END Genome()
        
        
        % STORE TESTED BRAIN %
        % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % %
        function insertBrain(G,brain)
            i = 1;
            while i <= size(G.newRank,1) && brain.fitness < G.newRank(i,2)
                i = i+1;
            end
            G.newRank = [G.newRank(1:i-1,:); brain.place, brain.fitness; G.newRank(i:end,:)];
            
            % Create new generation
            if size(G.newRank,1) == G.pop
                G.genCount = G.genCount + 1;
                G.curGen = G.newGen;
                G.rank = G.newRank;
                G.newRank = [];
                elite = ceil(G.elitism*G.pop);
                for j = 1:elite
                    n = size(G.rank,1)*(size(G.rank,1)-1)/2*rand();
                    p1 = 1;
                    t = G.pop;
                    while t < n
                        p1 = p1 + 1;
                        t = t + G.pop - p1 + 1;
                    end
                    G.newGen(j) = G.curGen(G.rank(p1,1));
                    G.newRank(j,:) = [j,G.rank(p1,2)];
                end
                G.nextBrain = elite + 1;
                save(strcat('genome_',datestr(datetime(),'yy_mm_dd_HH'),'.mat'),'G');
            end
        end
        
        % CREATES AND RETURNS A NEW BRAIN %
        % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % %
        function brain = newBrain(G)
            % pick 2 new genomes
            n = size(G.rank,1)*(size(G.rank,1)-1)/2*rand();
            p1 = 1;
            t = G.pop;
            while t < n
                p1 = p1 + 1;
                t = t + G.pop - p1;
            end
            n = size(G.rank,1)*(size(G.rank,1)-1)/2*rand();
            p2 = 1;
            t = G.pop;
            while t < n
                p2 = p2 + 1;
                t = t + G.pop - p2;
            end
            
            % Get actual index
            p1 = G.rank(p1,1);
            p2 = G.rank(p2,1);
    
            % CROSSOVER
            gene(G.height+1).m = 0;
            if G.height > 0
                rows = rand(G.width,1);
                gene(1).m = zeros(G.width,G.input+1);
                gene(1).m(rows<G.crossover,:)=G.curGen(p1).gene(1).m(rows<G.crossover,:);
                gene(1).m(rows>=G.crossover,:)=G.curGen(p2).gene(1).m(rows>=G.crossover,:);
                j = 2;
                while j <= G.height
                    rows = rand(G.width,1);
                    gene(j).m = zeros(G.width,G.width+1);
                    gene(j).m(rows<G.crossover,:)=G.curGen(p1).gene(j).m(rows<G.crossover,:);
                    gene(j).m(rows>=G.crossover,:)=G.curGen(p2).gene(j).m(rows>=G.crossover,:);
                    j = j+1;
                end
                rows = rand(G.output,1);
                gene(j).m = zeros(G.output,G.width+1);
                gene(j).m(rows<G.crossover,:)=G.curGen(p1).gene(j).m(rows<G.crossover,:);
                gene(j).m(rows>=G.crossover,:)=G.curGen(p2).gene(j).m(rows>=G.crossover,:);
            else
                rows = rand(G.output,1);
                gene(1).m = zeros(G.output,G.input+1);
                gene(1).m(rows<G.crossover,:)=G.curGen(p1).gene(1).m(rows<G.crossover,:);
                gene(1).m(rows>=G.crossover,:)=G.curGen(p2).gene(1).m(rows>=G.crossover,:);
            end
    
            % MUTABILITY
            for i = 1:G.height+1
                for j = 1:size(gene(i).m,1)
                    for k = 1:size(gene(i).m,2)
                        r = rand();
                        if r < G.mutability
                            gene(i).m(j,k) = 2*rand()-1;
                        elseif r < G.mutability*10
                            gene(i).m(j,k) = max(-1,min(1,gene(i).m(j,k)+(rand()-0.5)/10));
                        end
                    end
                end
            end
            
            n = G.nextBrain;
            G.nextBrain = G.nextBrain+1;
            G.newGen(n).gene = gene;
            G.newGen(n).fitness = 0;
            G.newGen(n).place = n;
        
            G.brainCount = G.brainCount + 1;
            brain = G.newGen(n);
        end % END NewBrain()
        
        function out = convergence(G)
            A = [];
            B = [];
            
            for a = 1:G.pop
                z = 1;
                gene = G.curGen(a).gene;
                for i = 1:G.height+1
                    for j = 1:size(gene(i).m,1)
                        for k = 1:size(gene(i).m,2);
                            B(z) = gene(i).m(j,k);
                            z = z + 1;
                        end
                    end
                end
                A = [A;B];
            end
            
            S = std(A);
            score = sum(S);
            percent = 100*(1 - score/(size(A,2)*sqrt(1/3)));
            out = [score,percent];
        end
    end
    methods (Static)
        
        % EVALUATE NET %
        function output = evaluateNet(net,input,p)
            output = input';
            for i = 1:size(net,2)
                output = sigmf(net(i).m*[output;-1],[p 0]);
            end       
        end % END evaluateNet()
        
    end % END METHODS
end

