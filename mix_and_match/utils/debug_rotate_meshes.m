function [Shapes_out] = debug_rotate_meshes(Shapes)

A = [0     0    -1
     0     1     0
     1     0     0];
 
Shapes_out = Shapes;
for i = 1:length(Shapes_out)
    Shapes_out{i}.vertexPoss = A*Shapes_out{i}.vertexPoss;
    Shapes_out{i} = normalize_shape(Shapes_out{i});
end

function [Shape] = normalize_shape(Shape_unnor)

%Normalize shape
pos_min = min(Shape_unnor.vertexPoss')';
pos_max = max(Shape_unnor.vertexPoss')';
scale = 1/(pos_max(2) - pos_min(2));

center = (pos_min+pos_max)/2;
numV = size(Shape_unnor.vertexPoss, 2);
Shape_unnor.vertexPoss = Shape_unnor.vertexPoss - center*ones(1, numV);
Shape_unnor.vertexPoss = Shape_unnor.vertexPoss*scale;

Shape = Shape_unnor;
 
