%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize the shape so that the center of the bounding box is at the
% origin, and the diagonal length of the bounding box is 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Shape] = mm_normalize_shape(Shape_unnor)

%Normalize shape
pos_min = min(Shape_unnor.vertexPoss')';
pos_max = max(Shape_unnor.vertexPoss')';
scale = 1/norm(pos_max - pos_min);

center = (pos_min+pos_max)/2;
numV = size(Shape_unnor.vertexPoss, 2);
Shape_unnor.vertexPoss = Shape_unnor.vertexPoss - center*ones(1, numV);
Shape_unnor.vertexPoss = Shape_unnor.vertexPoss*scale;

Shape = Shape_unnor;

nv = size(Shape.vertexPoss, 2);
nf = size(Shape.faceVIds, 2);

Shape.vertexPoss = single(zeros(3, nv));
Shape.faceVIds = uint32(zeros(3, nf));

off_v = 0;
off_f = 0;

for i = 1:length(Shape_unnor.meshes)
    mesh = Shape_unnor.meshes{i};
    nf_m = size(mesh.faceVIds, 2);
    nv_m = length(mesh.vertexIds);
    vertexPoss = Shape_unnor.vertexPoss(:, mesh.vertexIds);
    Shape.meshes{i}.vertexIds = (off_v+1):(off_v + nv_m);
    Shape.vertexPoss(:, Shape.meshes{i}.vertexIds) = vertexPoss;
    Shape.faceVIds(:, (off_f+1):(off_f+nf_m)) = mesh.faceVIds + off_v;
    off_v = off_v + nv_m;
    off_f = off_f + nf_m;
end