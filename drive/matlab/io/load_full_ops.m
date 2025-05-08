function [a0_full, b0_full, c0_full, u0_full, uk_full, mb, ms] =  load_full_ops(path)
% load_full_ops function loads the full ROM operators stored in the specified path
% specified by input variable path.
%
% Output:
% a0_full : Full stiffness matrix of size mb+1 x mb+1 (The +1 comes from the zeroth mode)
% b-1_full : Full   mass    matrix of size mb+1 x mb+1
% cu_full : Full advection tensor of size mb x mb+1 x mb+1
% u0_full : Vector of size mb+1, which contains the ROM coefficients of the
%           projection of initial conditons onto mb-reduced space
% uk_full : Matrix of size mb+1 x ns. Each column contains the ROM coefficients
%           of the projection of the snapshot onto mb-reduced space
% mb : total number of modes
% ns : nubmer of snapshots used to create those operators
% 
% This function can take times to load operators with nb >= 300.

   fprintf('Loading ROM operators and vectors... \n');
   fprintf('Currently only support velocity... \n');

   mb=dlmread(fullfile(path,"nb"));

   % load stiffness matrix
   a0_full = dlmread(fullfile(path,"au"));
   %size(a0_full)
   %mb+1
   %(mb+1)*(mb+1)
   a0_full = reshape(a0_full,mb+1,mb+1);

   % load mass matrix
   b0_full = dlmread(fullfile(path,"bu"));
   b0_full = reshape(b0_full,mb+1,mb+1);

   % load advection tensor
   c0_full = dlmread(fullfile(path,"cu"));
   c0_full = reshape(c0_full,mb,mb+1,mb+1);

   u0_full = dlmread(fullfile(path,"u0"));

   ms = dlmread(fullfile(path,"ns"));
   uk_full = dlmread(fullfile(path,"uk"));
   uk_full = reshape(uk_full,mb+1,ms);

   fprintf("done loading ... \n");
   
end

