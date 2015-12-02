%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initilaize the free-from deformation structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FFD] = sp_ffd_init_sym(Shape, gridRes)
% Initialize the FFD data structure, given a grid resolution
% Input arguments:
%   Shape: the input shape
%   gridRes: the resolution of the grid
% Ouput argument:
%   FFD: the initialize grid-structure

% Compute the bounding box
boxMin = double(min(Shape.vertexPoss'));
boxMax = double(max(Shape.vertexPoss'));
boxCenter = (boxMin + boxMax)/2;
boxSize = boxMax - boxMin;

% Initialize the grid
dimX = double(floor(boxSize(1)/gridRes) + 2);
dimY = double(floor(boxSize(2)/gridRes) + 2);
dimZ = double(floor(boxSize(3)/gridRes) + 2);

FFD.lc(1) = boxCenter(1) - gridRes*(dimX-1)/2;
FFD.lc(2) = boxCenter(2) - gridRes*(dimY-1)/2;
FFD.lc(3) = boxCenter(3) - gridRes*(dimZ-1)/2;
FFD.gridRes = gridRes;
FFD.dimX = dimX;
FFD.dimY = dimY;
FFD.dimZ = dimZ;

% Initilize the gridRes solution
FFD.ctrlPos_ori = zeros(3, dimX*dimY*dimZ);
IDX = 1:(dimX*dimY*dimZ);

offsetsY = floor((IDX-1)/(dimX*dimZ)) + 1;
offsetsZ = IDX - (offsetsY-1)*dimX*dimZ;
offsetsX = floor((offsetsZ-1)/dimZ) + 1;
offsetsZ = offsetsZ - (offsetsX-1)*dimZ;

FFD.ctrlPos_ori(1,:) = FFD.lc(1) + (offsetsX-1)*gridRes;
FFD.ctrlPos_ori(2,:) = FFD.lc(2) + (offsetsY-1)*gridRes;
FFD.ctrlPos_ori(3,:) = FFD.lc(3) + (offsetsZ-1)*gridRes;
FFD.ctrlPos_cur = FFD.ctrlPos_ori;

FFD.offsetsX = offsetsX;
FFD.offsetsY = offsetsY;
FFD.offsetsZ = offsetsZ;

if mod(dimX, 2) == 0
    colids_1 = find(offsetsX < (dimX+1)/2);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
        (dimX - offsetsX(colids_1))*dimZ +...
        offsetsZ(colids_1);
    
    rows = kron(1:length(colids_1), ones(1,2));
    cols = reshape([colids_1;colids_2], [1, 2*length(colids_1)]);
    vals = kron(ones(1,length(colids_1)), [1,1]);
    FFD.Ax = sparse(rows, cols, vals);
    FFD.bx = zeros(size(FFD.Ax,1),1);
    
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    FFD.Az = sparse(rows, cols, vals);
    FFD.bz = zeros(size(FFD.Az,1),1);
    
    colids_1 = find(offsetsY == 1);
    colids_2 = find(offsetsY == dimY);
    rows = 1:(2*dimX*dimZ);
    cols = double([colids_1, colids_2]);
    vals = ones(1, 2*dimX*dimZ);
    Ay1 = sparse(rows, cols, vals, 2*dimX*dimZ, dimX*dimY*dimZ);
    yMin = FFD.lc(2);
    yMax = yMin + (dimY-1)*gridRes;
    by1 = [ones(1, length(colids_1))*yMin, ones(1, length(colids_2))*yMax]';
    
    colids_1 = find(offsetsX < (dimX+1)/2 & offsetsY > 1 & offsetsY < dimY);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
    (dimX - offsetsX(colids_1))*dimZ +...
    offsetsZ(colids_1);
    rows = kron(1:length(colids_1), ones(1,2));
    cols = reshape([colids_1;colids_2], [1, 2*length(colids_1)]);
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    Ay2 = sparse(rows, cols, vals, length(colids_1), dimX*dimY*dimZ);
    by2 = zeros(size(Ay2,1),1);
    
    FFD.Ay = [Ay1;Ay2];
    FFD.by = [by1;by2];
else
    cols = find(offsetsX == (dimX+1)/2);
    rows = 1:length(cols);
    vals = ones(1, length(cols));
    
    colids_1 = find(offsetsX < dimX/2);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
        (dimX - offsetsX(colids_1))*dimZ +...
        offsetsZ(colids_1);
    
    rows = [rows, length(rows)+kron(1:length(colids_1), ones(1,2))];
    cols = [cols, reshape([colids_1;colids_2], [1, 2*length(colids_1)])];
    vals = [vals, kron(ones(1,length(colids_1)), [1,1])];
    FFD.Ax = sparse(rows, cols, vals);
    FFD.bx = zeros(size(FFD.Ax,1),1);
    
    rows = [kron(1:length(colids_1), ones(1,2))];
    cols = [reshape([colids_1;colids_2], [1, 2*length(colids_1)])];
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    FFD.Az = sparse(rows, cols, vals);
    FFD.bz = zeros(size(FFD.Az,1),1);
    
    colids_1 = find(offsetsY == 1);
    colids_2 = find(offsetsY == dimY);
    rows = 1:(2*dimX*dimZ);
    cols = double([colids_1, colids_2]);
    vals = ones(1, 2*dimX*dimZ);
    Ay1 = sparse(rows, cols, vals, 2*dimX*dimZ, dimX*dimY*dimZ);
    yMin = FFD.lc(2);
    yMax = yMin + (dimY-1)*gridRes;
    by1 = [ones(1, length(colids_1))*yMin, ones(1, length(colids_2))*yMax]';
    
    colids_1 = find(offsetsX < dimX/2 & offsetsY > 1 & offsetsY < dimY);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
    (dimX - offsetsX(colids_1))*dimZ +...
    offsetsZ(colids_1);
    rows = kron(1:length(colids_1), ones(1,2));
    cols = reshape([colids_1;colids_2], [1, 2*length(colids_1)]);
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    Ay2 = sparse(rows, cols, vals, length(colids_1), dimX*dimY*dimZ);
    by2 = zeros(size(Ay2,1),1);
    
    FFD.Ay = [Ay1;Ay2];
    FFD.by = [by1;by2];
end

% The smoothness term
subset = find(offsetsX < dimX -1);
vIds0 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds1 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset))*dimZ + offsetsZ(subset);
vIds2 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)+1)*dimZ + offsetsZ(subset);

rows = kron(1:length(vIds0), ones(1,3));
cols = [vIds0;vIds1;vIds2];
vals = kron(ones(1,length(vIds0)), [1,-2,1]);
Jx = sparse(rows, cols, vals, length(vIds0), dimX*dimY*dimZ);

subset = find(offsetsY < dimY -1);
vIds0 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds1 = (offsetsY(subset))*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds2 = (offsetsY(subset)+1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);

rows = kron(1:length(vIds0), ones(1,3));
cols = [vIds0;vIds1;vIds2];
vals = kron(ones(1,length(vIds0)), [1,-2,1]);
Jy = sparse(rows, cols, vals, length(vIds0), dimX*dimY*dimZ);

subset = find(offsetsZ < dimZ -1);
vIds0 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds1 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset)+1;
vIds2 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset)+2;

rows = kron(1:length(vIds0), ones(1,3));
cols = [vIds0;vIds1;vIds2];
vals = kron(ones(1,length(vIds0)), [1,-2,1]);
Jz = sparse(rows, cols, vals, length(vIds0), dimX*dimY*dimZ);

J_smooth = [Jx;Jy;Jz];
FFD.H_smooth = J_smooth'*J_smooth;
