function [shapeDefmScore, verDefmScores] = i2s_eval_deformation(...
    Shape, FFD, numSamples, knn)
% Compute the local and global distortions of a FFD
% Input arguments:
%   Shape: the input shape
%   FFD: the free-form deformation to be evaluated
%   numSamples: as vertices of Warehouse shapes are irregularly sampled, we
%               augment the point cloud (the vertex positions) with uniform
%               samples on the point cloud  
%   knn: the neighborhood size used in evaluating local distortiion scores
% Output arguments:
%   shapeDefmScore: the global deformation score
%   verDefmScores: the local deformation scores

% Compute deform vertices and original vertices
vertexPoss_ori = double(Shape.vertexPoss);
points = i2s_mesh_sampling(Shape, numSamples);
% only keep the points who distances to vertexPoss_ori is less than a
% threshold
[IDX, dis] = knnsearch(vertexPoss_ori', points');
ids = find(dis > 1e-2);
points_augmented = [vertexPoss_ori, points(:, ids)];
%
ffd_coeff_vertexPoss = i2s_ffd_basis_coeff(FFD, points_augmented);
points_augmented_deform = FFD.ctrlPos_cur*ffd_coeff_vertexPoss';
%
shapeDefmScore = eval_distortion(points_augmented,...
    points_augmented_deform);

numV = size(Shape.vertexPoss, 2);
IDX = knnsearch(points_augmented',points_augmented(:, 1:numV)', 'k', knn);

verDefmScores = zeros(1, numV);
for vId = 1 : numV
    ids = IDX(vId, :);
    verDefmScores(vId) = eval_distortion(points_augmented(:, ids),...
        points_augmented_deform(:,ids));
end


function [distortion] = eval_distortion(points_ori, points_deformed)
%
numP = size(points_ori, 2);
[R, t] = horn87(points_ori, points_deformed, ones(1, numP));
d = R*points_ori + t*ones(1, numP) - points_deformed;
distortion = sqrt(sum(sum(d.*d))/numP);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform uniform sampling of each mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points] = i2s_mesh_sampling(Shape, numSamples)
%
numFaces = size(Shape.faceVIds, 2);
faceAreas = zeros(1, numFaces);
for faceId = 1:numFaces
    p1 = Shape.vertexPoss(:, Shape.faceVIds(1, faceId));
    p2 = Shape.vertexPoss(:, Shape.faceVIds(2, faceId));
    p3 = Shape.vertexPoss(:, Shape.faceVIds(3, faceId));
    e12 = p1 - p2;
    e13 = p1 - p3;
    e23 = p2 - p3;
    a = e12'*e12;
    b = e13'*e13;
    c = e23'*e23;
    faceAreas(faceId) = sqrt(2*(a*b+a*c+b*c) - (a*a+b*b+c*c))/4;
    if faceId > 1
        faceAreas(faceId) = faceAreas(faceId) + faceAreas(faceId-1);
    end
end

faceAreas = faceAreas/faceAreas(numFaces);
sample_ts = sort(rand(1, numSamples));
points = zeros(3, numSamples);

faceId = 1;
for sId = 1 : numSamples
    while sample_ts(sId) > faceAreas(faceId)
        faceId = faceId + 1;
    end
    p1 = Shape.vertexPoss(:, Shape.faceVIds(1,faceId));
    p2 = Shape.vertexPoss(:, Shape.faceVIds(2,faceId));
    p3 = Shape.vertexPoss(:, Shape.faceVIds(3,faceId));
    r1 = rand(1,1);
    r2 = rand(1,1);
    t1 = (1-sqrt(r1));
    t2 = sqrt(r1)*(1-r2);
    t3 = sqrt(r1)*r2;
    points(:, sId) = t1*p1 + t2*p2 + t3*p3;
end


