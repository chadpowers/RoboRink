classdef TrainBrain < handle
    %GENOME, this class contains the entire gene pool for the neural net
    % AI that controls the computer mallet, genome can be shared across
    % various, simultaneous simulations of mallet
    
    properties (SetAccess = private)
        % NETWORK PROPERTIES
        height;       % height of hidden layers
        width;        % width of hidden layers
        input;        % number of inputs
        output;       % number of outputs
        learningrate; % learning rate
        
        % EVOLUTION PROPERTIES
        
        % HELPER/COUNTERS
        genCount; % number of generations
        brainCount;      % number of brains given
        writefile;  %file to write too
        
        convergenceTrend; % for plot
        topRawFitnessTrend; % for plot
        medRawFitnessTrend;
        botRawFitnessTrend;
        
        brain;  % 

    end
    
    properties
        mutability; % percent chance of complete mutation - 10x for minor mutation
    end
    
    properties (Access = private)
        nextBrain;
        newGen;  % all genes of old generation
    end
    
    methods
        
        % CONSTRUCTOR %
        function G = TrainBrain(h,w,in,out,learn,wf)
            % INITIALIZE VARIABLES
            G.height        = h;
            G.width         = w;
            G.input         = in;
            G.output        = out;
            G.learningrate  = learn;
            G.writefile     = wf;
            
            G.convergenceTrend = []; % for plot
            G.topRawFitnessTrend = []; % for plot
            G.medRawFitnessTrend = [];
            G.botRawFitnessTrend = [];
            
            if G.height > 0
                gene(1).m = 2*(rand(G.width,G.input+1)-0.5);
                j = 2;
                while j <= G.height
                    gene(j).m = 2*(rand(G.width,G.width+1)-0.5);
                    j = j+1;
                end
                gene(j).m = 2*(rand(G.output,G.width+1)-0.5);
            else
                gene(1).m = 2*(rand(G.output,G.input+1)-0.5);
            end
            G.brain.gene = gene;
            
        end % END Genome()
        
        function learn(G,input,d_output)
            [b_inputs,b_outputs] = G.evaluateNet(G.brain.gene,input,1);
            d(G.height+1).v=d_output-b_outputs(G.height+1).v;
            for i = 0:G.height
                j = G.height +1 - i;
                err(j).v = (d(j).v).*b_outputs(j).v.*(1-b_outputs(j).v);
                G.brain.gene(j).m = G.brain.gene(j).m+G.learningrate*err(j).v*b_inputs(j).v;
                if j > 1
                    d(j-1).v = G.brain.gene(j).m(:,1:end-1)'*err(j).v;
                end
            end
        end
        
        
    end
    methods (Static)
        
        % EVALUATE NET %
        function [inputs,outputs] = evaluateNet(net,input,p)
            inputs(1).v = [input,-1];
            for i = 1:size(net,2)
                outputs(i).v = sigmf(net(i).m*inputs(i).v',[p 0]);
                if i+1 <= size(net,2)
                    inputs(i+1).v = [outputs(i).v',-1];
                end
            end       
        end % END evaluateNet()
        
    end % END METHODS
end

