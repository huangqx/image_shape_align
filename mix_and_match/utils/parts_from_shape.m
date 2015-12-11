function [meshes] = mm_parts_from_shape(Shape, eps)
%
numVertices = size(Shape.vertexPoss, 2);
ids = zeros(1, numVertices);

id = 1;
for i = 1:numVertices
    if ids(i) > 0
        continue;
    end
    d = Shape.vertexPoss(:, i)*single(ones(1, numVertices)) - Shape.vertexPoss;
    tp = find(sqrt(sum(d.*d)) < eps);
    tp = tp(find(tp > i));
end