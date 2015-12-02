%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform uniform sampling of each mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points] = sp_mesh_sampling(Shape, numSamples)
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
