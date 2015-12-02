%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the coefficient of 3D positions with respect to a free-form
% deformation structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = sp_ffd_basis_coeff(FFD, poss)
%
poss = double(poss);

tX = (poss(1,:) - FFD.lc(1))/FFD.gridRes;
idX = floor(tX) + 1;
tX = tX - (idX-1);

tY = (poss(2,:) - FFD.lc(2))/FFD.gridRes;
idY = floor(tY) + 1;
tY = tY - (idY-1);

tZ = (poss(3,:) - FFD.lc(3))/FFD.gridRes;
idZ = floor(tZ) + 1;
tZ = tZ - (idZ-1);

nP = size(poss, 2);

rows = ones(8,1)*(1:nP);
cols = zeros(8, nP);
vals = zeros(8, nP);

cols(1,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ;
vals(1,:) = (1-tX).*(1-tY).*(1-tZ);

cols(2,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ + 1;
vals(2,:) = (1-tX).*(1-tY).*tZ;

cols(3,:) = (idY)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ;
vals(3,:) = (1-tX).*tY.*(1-tZ);

cols(4,:) = (idY)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ + 1;
vals(4,:) = (1-tX).*tY.*tZ;

cols(5,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ;
vals(5,:) = tX.*(1-tY).*(1-tZ);

cols(6,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ + 1;
vals(6,:) = tX.*(1-tY).*tZ;

cols(7,:) = (idY)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ;
vals(7,:) =  tX.*tY.*(1-tZ);

cols(8,:) = (idY)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ + 1;
vals(8,:) =  tX.*tY.*tZ;

C = sparse(rows, cols, vals, nP, FFD.dimX*FFD.dimY*FFD.dimZ);