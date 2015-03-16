function [z e] = SimakovFrolovaBasri2003(Io, mo, Ij, Mj, Tj, p)
%
% Re-implementation of 
% Simakov, Frolova and Basri "Dense Shape Reconstruction from a Moving
% Oject under Arbitrary, Unknown Lighting", 2003
%
%
% This re-imp is restricted to linear lighting model (eq. 10 in the  paper)
%
%
% Usage:
%   z = SimakovFrolovaBasri2003(Io, mo, Ij, Mj, Tj, p)
%
% Inputs:
%   Io  - reference image, 2D matrix of intensity values
%   mo  - mask for reference image, 2D binary image of the same size as Io.
%         Depth is estimated only for pixels in the mask.
%   Ij  - cell array of corresponding images (2D matrices)
%   Mj  - cell array of corresponding masks (2D matrices) same size as Ij's
%   Tj  - struct array of transformations, each element has
%         ().R  - rotation matrix w.r.t reference frame (3x3 matrix)
%         ().t  - translation vector w.r.t reference frame (2x1 vector)
%         ().scale - scaling scalar w.r.t reference frame 
%   p   - parameters structure
%         ().SFB_or_BC - using SFB consistency measure, or
%                        brightness constancy (BC)
%         ().z_range   - [min_z max_z] range of labels
%         ().NL        - number of labels (100 - 500)
%         ().L1_truncate - percent of Z range (.25 - .5)
%         ().lambda    - smoothness cost weight (roughly 1e-4 for 'bc', 5e-3 for 'sfb')
%
% Output:
%   z   - depth values for pixels in Io
%   e   - energy of solution
%
%
%
% 
% 
%
% Copyright (c) Bagon Shai
% Department of Computer Science and Applied Mathmatics
% Wiezmann Institute of Science
% http://www.wisdom.weizmann.ac.il/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
%
%
% If used in an academic research, the follwoing citation must be included
% in any resulting publication:
%
%   [1] Denis Simakov, Darya Frolova and Ronen Basri 
%       "Dense Shape Reconstruction from a Moving Oject under Arbitrary, Unknown Lighting", 
%       ICCV, 2003
% 
%   [2] Shai Bagon 
%       "Matlab Implementation of Simakov Frolova and Basri 3D Reconstruction",
%       June, 2012
%
%
% June. 2012
%
% 

INF = 1e6;

%------------
% determine the label space
Z = linspace(p.z_range(1), p.z_range(2), p.NL);
num_of_labels = numel(Z);


num_of_nodes = nnz(mo);
%
mapping_image_nodes = zeros(numel(mo),1);
mapping_image_nodes(mo) = 1:num_of_nodes;
mapping_nodes_image = find(mapping_image_nodes > 0);


%-------------
% Computing data term according to SFB or BC
DataCost = zeros(num_of_nodes,num_of_labels);


% estimate the matrix Omega (all rotations)
NJ = numel(Ij);

switch upper(p.SFB_or_BC)
    case {'SFB'}
        Omega = zeros(NJ,3);
        for ji=1:NJ
            [theta w] = AxisAngleRep(Tj(ji).R);
            Omega(ji,:) = -theta * w';
        end
        piOmega = Omega * inv(Omega'*Omega) * Omega' - 1;
    case {'BC'}
        piOmega = eye(NJ);
    otherwise
        error('SimakovFrolovaBasri2003:cost','unknown cost %s', p.SFB_or_BC);
end

% define grid size
sz = size(mo); %sz1 denotes #rows, sz2 denotes #cols
sz1=sz(1);
sz2=sz(2);
[xexpand,yexpand] = meshgrid(1:sz2,1:sz1);


% for each depth value Z(li)
for li=1:num_of_labels
    fprintf(1, 'Dc for depth %.1f (%d/%d)\n', Z(li), li, num_of_labels);
    
    dI = zeros(NJ,num_of_nodes);
    
    % for each J image
    for ji = 1:NJ
        
        Rot_scale = Tj(ji).R * Tj(ji).scale;
        Trans = Tj(ji).t;
        

        xtag = xexpand*Rot_scale(1,1) + yexpand*Rot_scale(1,2) + ...
            Z(li)*Rot_scale(1,3) + Trans(1);
        ytag = xexpand*Rot_scale(2,1) + yexpand*Rot_scale(2,2) + ...
            Z(li)*Rot_scale(2,3) + Trans(2);

        J = Ij{ji};
        J(~Mj{ji}) = INF;
        Jinterp = interp2(xexpand,yexpand,J,xtag,ytag,'*linear');
        Jinterp(isnan(Jinterp)) = INF;
        
        
        dI(ji,:) = Io(mo)-Jinterp(mo);
        
    end % for each J
    
    muP = (piOmega * dI)';
    
    DataCost(:,li) = sum(muP.*muP, 2);
end % for each depth Z
        

%--------
% optimization part
fprintf(1, 'optimizing...\n');

% grid
sz = size(Io);
[ii_image jj_image] = sparse_adj_matrix(sz, 1, 1, 1);
% sel = jj>ii; % the right direction... (both x and y)
% ii_image=ii(sel);
% jj_image= jj(sel);
ii_nodes = mapping_image_nodes(ii_image);
jj_nodes = mapping_image_nodes(jj_image);
index_nodes = (ii_nodes > 0) & (jj_nodes > 0);
ii_nodes = ii_nodes(index_nodes);
jj_nodes = jj_nodes(index_nodes);

Sparse_Graph = sparse(ii_nodes, jj_nodes, 1, num_of_nodes, num_of_nodes);

% truncated L1 smoothness cost
V = abs( bsxfun(@minus, Z, Z') );
V = p.lambda * min( V, p.L1_truncate * abs(Z(end)-Z(1)) );

[Ubias, wta] = min(DataCost, [], 2); % Winners Takes All inititialization
DataCost = bsxfun(@minus, DataCost, Ubias);


gch = GraphCut('open', DataCost', V, Sparse_Graph);
gch = GraphCut('t',gch,true); % truncate due to numerical precissions
gch = GraphCut('set', gch, wta-1);
[gch x] = GraphCut('expand', gch);
x = double(x)+1;
[gch e] = GraphCut('energy',gch);
gch = GraphCut('close', gch);

z = zeros(sz);
z(mapping_nodes_image) = Z(x);

