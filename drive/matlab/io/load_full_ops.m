function [a0_full, b0_full, c0_full, u0_full, uk_full, mb, ms] =  load_full_ops(path)
% load_full_ops loads the full ROM operators stored in path.
%
% Output:
% - a0_full: full stiffness matrix of size mb+1 x mb+1. The +1 comes from the zeroth mode.
% - b0_full: full mass matrix of size mb+1 x mb+1.
% - c0_full: full advection tensor of size mb x mb+1 x mb+1.
% - u0_full: vector of size mb+1 containing the ROM coefficients of the projection of the initial condition.
% - uk_full: matrix of size mb+1 x ns. Each column contains the ROM coefficients of the projection of one snapshot.
% - mb: total number of modes.
% - ms: number of snapshots used to create the operators.
%
% This function can take time to load operators when nb >= 300.

   fprintf('Loading ROM operators and vectors... \n');
   fprintf('Currently only support velocity... \n');

   mb=dlmread(fullfile(path,'nb'));

   % load stiffness matrix
   a0_full = dlmread(fullfile(path,'au'));
   %size(a0_full)
   %mb+1
   %(mb+1)*(mb+1)
   a0_full = reshape(a0_full,mb+1,mb+1);

   % load mass matrix
   b0_full = dlmread(fullfile(path,'bu'));
   b0_full = reshape(b0_full,mb+1,mb+1);

   % load advection tensor
   c0_full = dlmread(fullfile(path,'cu'));
   c0_full = reshape(c0_full,mb,mb+1,mb+1);

   u0_full = dlmread(fullfile(path,'u0'));

   ms = dlmread(fullfile(path,'ns'));
   uk_full = dlmread(fullfile(path,'uk'));
   uk_full = reshape(uk_full,mb+1,ms);

   fprintf('done loading ... \n');
   
end
