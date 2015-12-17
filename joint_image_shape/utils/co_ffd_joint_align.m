function [patch_ffds] = co_ffd_joint_align(patches, maps, Para)
% This functions optimizes a free-from deformation for each image patch
% to align them in a canonical domain (the one given by the first image
% patch). Note that we use reweighted least squares to optimize the L1-norm
% Input arguments:
%  patches: a cell array of images
%  maps: the pre-computed dense correspondences between pairs of images
%  Para: Para.gridRes_2d_ffd_alignment: resolution of the 2D ffd domain
%        Para.lambda_2d_ffd_alignment: the smoothness weight of the
%        alignment
% Output argument:
%  patch_ffds: the optimized free-form deformations


numPatches = length(patches);
offsets = zeros(1, numPatches+1);
% used in accessing ctrl points of each patch

Basis = cell(1, numPatches);
Hs = cell(1, numPatches);
ctrlPts = cell(1, numPatches);
for i = 1:numPatches
    [Basis{i}, Hs{i}, ctrlPts{i}] = ffd_basis(patches{i},...
        Para.gridRes_2d_ffd_alignment);
    offsets(i+1) = offsets(i) + size(Basis{i}, 2);
end

% store the original ctrl points
x_ori = zeros(2, offsets(numPatches+1));
for i = 1:numPatches
    ids = (offsets(i)+1):offsets(i+1);
    x_ori(:, ids) = ctrlPts{i};
end

numCorresAll = 0;
for id = 1:length(maps)
    % Allocate space for correspondence weights which will be used to 
    % reweighted least squares
    maps{id}.w_st = ones(1, size(maps{id}.corres_st, 2));
    maps{id}.w_ts = ones(1, size(maps{id}.corres_ts, 2));
    numCorresAll = numCorresAll + length(maps{id}.w_st);
    numCorresAll = numCorresAll + length(maps{id}.w_ts);
end

nCtrls = offsets(length(offsets));

for outerIter = 1:16
    % Optimize the control points
    A = sparse(nCtrls, nCtrls);
    
    % Data term
    for id = 1:length(maps)
        map = maps{id};
        rIds = [(offsets(map.sId)+1):offsets(map.sId+1),...
                (offsets(map.tId)+1):offsets(map.tId+1)];
        Jst = [Basis{map.sId}(map.corres_st(1,:)',:),...
              -Basis{map.tId}(map.corres_st(2,:)',:)];
        Wst = sparse(1:length(map.w_st), 1:length(map.w_st), map.w_st);
        Jts = [Basis{map.sId}(map.corres_ts(2,:)',:),...
              -Basis{map.tId}(map.corres_ts(1,:)',:)];
        Wts = sparse(1:length(map.w_ts), 1:length(map.w_ts), map.w_ts);
        B = Jst'*Wst*Jst + Jts'*Wts*Jts;
        A(rIds, rIds) = A(rIds, rIds) + B*map.weight;
    end
    
    % regularization term
    for id = 1:numPatches
        rIds = (offsets(id)+1):offsets(id+1);
        A(rIds, rIds) = A(rIds, rIds) + Para.lambda_2d_ffd_alignment*Hs{id};
    end
    
    % Manipulate ctrl points of the first patch
    fixed_ctrls = ctrlPts{1};
    ids_0 = 1:offsets(2);
    ids_1 = (offsets(2)+1):nCtrls;
    A11 = A(ids_1, ids_1);
    b = - A(ids_1, ids_0)*fixed_ctrls';
    
    n = length(ids_1);
    A11 =  A11 + Para.lambda_2d_ffd_alignment*sparse(1:n, 1:n, ones(1,n))/8;
    b = b + Para.lambda_2d_ffd_alignment*x_ori(:, ids_1)'/8;
    
    remaining_ctrls = A11\b;
    for i = 2:length(ctrlPts)
        rIds = ((offsets(i)+1):offsets(i+1)) - offsets(2);
        ctrlPts{i} = remaining_ctrls(rIds,:)';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reweight the initial correspondences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tp = zeros(1, numCorresAll);
    off = 0;
    for id = 1:length(maps)
        map = maps{id};
        dst = Basis{map.sId}(map.corres_st(1,:)',:)*ctrlPts{map.sId}'...
            - Basis{map.tId}(map.corres_st(2,:)',:)*ctrlPts{map.tId}';
        maps{id}.w_st = sum(dst'.*dst');
        dts = Basis{map.sId}(map.corres_ts(2,:)',:)*ctrlPts{map.sId}'...
            - Basis{map.tId}(map.corres_ts(1,:)',:)*ctrlPts{map.tId}';
        maps{id}.w_ts = sum(dts'.*dts');
        tp((off+1):(off+length(maps{id}.w_st))) = maps{id}.w_st;
        off = off + length(maps{id}.w_st);
        tp((off+1):(off+length(maps{id}.w_ts))) = maps{id}.w_ts;
        off = off + length(maps{id}.w_ts);
    end
    sigma2 = median(tp);
    for id = 1:length(maps)
        maps{id}.w_st = sigma2./(maps{id}.w_st+sigma2);
        maps{id}.w_ts = sigma2./(maps{id}.w_ts+sigma2);
    end
end

patch_ffds = cell(1, numPatches);
for id = 1:numPatches
    patch_ffds{id}.ctrlPts = ctrlPts{id};
    patch_ffds{id}.basis = Basis{id};
end


function [B, Hs, ctrlPts] = ffd_basis(patch, gridRes)

[nRows, nCols, k] = size(patch);
cellWidth = nCols/gridRes;
cellHeight = nRows/gridRes;
B = zeros(nRows*nCols, (gridRes+1)*(gridRes+1));

ctrlPts = [kron(0:gridRes, ones(1,gridRes+1))*cellWidth;
    kron(ones(1,gridRes+1),0:gridRes)*cellHeight];

for col = 1 : nCols
    for row = 1 : nRows
        pointId = row + (col-1)*nRows;
        x1 = (col-0.5)/cellWidth;
        y1 = (row-0.5)/cellHeight;
        
        ctrlP_colId = floor(x1) + 1;
        ctrlP_rowId = floor(y1) + 1;
        x = x1 + 1 - ctrlP_colId;
        y = y1 + 1 - ctrlP_rowId;
        
        ctrlP_id00 = (ctrlP_colId-1)*(gridRes+1) + ctrlP_rowId;
        ctrlP_id01 = ctrlP_colId*(gridRes+1) + ctrlP_rowId;
        ctrlP_id10 = (ctrlP_colId-1)*(gridRes+1) + ctrlP_rowId + 1;
        ctrlP_id11 = ctrlP_colId*(gridRes+1) + ctrlP_rowId + 1;
        B(pointId, ctrlP_id00) = (1-x)*(1-y);
        B(pointId, ctrlP_id01) = x*(1-y);
        B(pointId, ctrlP_id10) = (1-x)*y;
        B(pointId, ctrlP_id11) = x*y;
    end
end

J = zeros(2*(gridRes-1)*(gridRes+1), (gridRes+1)*(gridRes+1));
off = 1;
for cId = 1:(gridRes+1)
    for rId = 1:(gridRes-1)
        id0 = (cId-1)*(gridRes+1) + rId;
        id1 = (cId-1)*(gridRes+1) + rId+1;
        id2 = (cId-1)*(gridRes+1) + rId+2;
        J(off, id0) = 1;
        J(off, id1) = -2;
        J(off, id2) = 1;
        off = off + 1;
    end
end

for cId = 1:(gridRes-1)
    for rId = 1:(gridRes+1)
        id0 = (cId-1)*(gridRes+1) + rId;
        id1 = (cId+0)*(gridRes+1) + rId;
        id2 = (cId+1)*(gridRes+1) + rId;
        J(off, id0) = 1;
        J(off, id1) = -2;
        J(off, id2) = 1;
        off = off + 1;
    end
end

Hs = J'*J;