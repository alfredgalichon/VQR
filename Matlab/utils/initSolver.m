function [ solverR ] = initSolver( solver,path )
%INITSOLVER Summary of this function goes here
%   Detailed explanation goes here
if solver(1)=='g' | solver(1)=='G'
    solverR='gurobi';
    currentFolder = pwd;
    cd(path);
    gurobi_setup;
    cd(currentFolder);
else if solver(1)=='k' | solver(1)=='K'
        solverR='knitro';
    end
end

end

