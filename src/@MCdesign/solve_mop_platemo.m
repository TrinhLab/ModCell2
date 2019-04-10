function [x,fval,output,population,scores] = solve_mop_platemo(obj, options)
% Solves MOP with PlatEMO MOEAs.  Input and output mimics Matlab function gamultiobj()
% 
% Notes:
%   *If memory usage is high, note that platemo stores intermediate solutions
%   * The GLOBAL.m class in platemo was modified to include mcdesign as a property.

%The modcell.m problem in platemo needs the information below:
obj.ga_parameters.initial_population    = options.InitialPopulationMatrix;
obj.ga_parameters.n_objectives          = obj.prodnet.n_prod;
obj.ga_parameters.n_variables           = size(options.InitialPopulationMatrix,2);

n_evaluations = obj.ga_parameters.stall_generations*obj.ga_parameters.population_size;
output.generations = obj.ga_parameters.stall_generations;

%% Platemo
Global          = GLOBAL('-algorithm',str2func(obj.ga_parameters.algorithm),...
                    '-problem',@modcell, '-evaluation',n_evaluations,'-mcdesign',obj,...
                    '-mode', 3, '-outputFcn', @(obj)[]);
Global.mcdesign = obj;
%solve
Global.Start(); 

%% Retrieve raw population as well as pareto set and pareto front
warning off MATLAB:structOnObject
g           = struct(Global); % hack to access private properties
Population  = g.result{end}; % last Population

population  = zeros(length(Population),obj.ga_parameters.n_variables);
scores      = zeros(length(Population),obj.ga_parameters.n_objectives);

for i =1:length(Population)
    population(i,:) = Population(i).dec;
    scores(i,:)     = Population(i).obj;
end

%Get points in first non-dominated front
FrontNo = NDSort(scores,1);
x       = population(FrontNo==1,:);
fval    = scores(FrontNo==1,:);
