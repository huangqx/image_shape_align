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