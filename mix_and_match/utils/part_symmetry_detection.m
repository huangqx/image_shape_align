function [sym_pairs, nonsym_ids] = part_symmetry_detection(Shape)
% This function detects the symmetry relations among parts (disconnected
% components)
% 
Shape = mm_normalize_shape(Shape);

% Compute the coordinate system of each part
numParts = length(Shape.meshes);
for i = 1:length(Shape.meshes)
    part_frames{i} = first_order_moment(Shape, Shape.meshes{i});
end

%
nonsym_ids = [];
sym_pairs = [];
flags = ones(1, numParts);
for  i = 1 : numParts
    if flags(i) == 0
        continue;
    end
    minDif = 10;
    reflect_partId = -1;
    for j = (i+1):numParts
        if flags(j) == 0
            continue;
        end
        sizeDif = part_frames{i}.pca - part_frames{j}.pca;
        posDif = (part_frames{i}.center - part_frames{j}.center);
        posDif(1) = (part_frames{i}.center(1) + part_frames{j}.center(1))/2;
        dif = sqrt(sizeDif'*sizeDif + posDif'*posDif);
        if dif < minDif
            minDif = dif;
            reflect_partId = j;
        end
    end
    if minDif < 3e-2 && abs(part_frames{i}.center(1)) > 2e-2
        flags(reflect_partId) = 0;
        flags(i) = 0;
        sym_pairs = [sym_pairs, [i; reflect_partId]];
    else
        nonsym_ids = [nonsym_ids, i];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the first-order moment of a part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame] = first_order_moment(Shape, mesh)

vertexPoss = Shape.vertexPoss(:, mesh.vertexIds);

m0 = 0;
m1 = zeros(3,1);
M2 = zeros(3,3);

for fId = 1 : size(mesh.faceVIds, 2)
    v1 = vertexPoss(:, mesh.faceVIds(1, fId));
    v2 = vertexPoss(:, mesh.faceVIds(2, fId));
    v3 = vertexPoss(:, mesh.faceVIds(3, fId));
    a = norm(cross(v2-v1, v3-v1));
    V = [v1, v2, v3]';
    S = [2,1,1;1,2,1;1,1,2]/24;
    C = a*V'*S*V;
    M2 = M2 + C;
    m1 = m1 + (v1+v2+v3)*a/6;
    m0 = m0 + a/2;
end

frame.center = m1/m0;

M2 = M2/m0;
m1 = m1/m0;
M2 = M2 - m1*m1';

[U,V] = eig(M2);
[s,ids] = sort(diag(V));
U = U(:, ids);
V = V(ids, ids);
frame.pca = sqrt(diag(V));
frame.axis = U;
frame.center = m1;

