%
% Main for running the simulation as a Matlab function.
%

function [] = main_matlab( parameter_file, runID )

if (nargin < 2)
    fprintf('Error: Missing required arguments: parameters file, run ID.\n');
    exit
end    
addpath( genpath('./FEM/Matlab') );     % low-level FEM implementation
run( parameter_file );                  % load parameters struct P
P.runID = runID;

%
% Initialize results folder.
%
dt = datetime;
results_folder = sprintf('../%s_%d-%.2d-%.2d--%.2d-%.2d-%.2d', ... 
                         P.runID, dt.Year, dt.Month, dt.Day, ... 
                         dt.Hour, dt.Minute, floor(dt.Second));
if (exist(results_folder, 'dir') == 7)
    % Folder exists; append process ID to the name.
    results_folder = join([results_folder, '_', num2str(pid)]);
end
[~,~,mID] = mkdir(results_folder);
if (~isempty(mID))
    fprintf('Error: Failed to create folder ''%s''.\n', results_folder);
    exit
end

copyfile(parameter_file, results_folder)

%
% Enter simulation.
%
Algorithm1( P, results_folder );

end     % END main_matlab().
