function [SymPart] = part_symmetry_detection(Shape)
% This function detects the symmetry relations among parts (disconnected
% components)
% 
Shape = normalize_shape(Shape);

%
numParts = max(Shape.meshes);

off = 0;
for partId = 1 : numParts
    faceIds = find(shape.partIds == partId);
    if length(faceIds) >= 2
        off = off + 1;
        [mesh, frame] = extract_shape_part(shape, faceIds);
        parts{off}.mesh = mesh;
        parts{off}.frame = frame;
        posY(off) = frame.center(2);
    end
end
[s, ids] = sort(posY);
parts = parts(ids);
parts1 = parts;

numParts = length(parts);
flags = ones(1, numParts);
parts_pairs = [];
parts_non_pairs = [];

for partId = 1 : numParts
    if partId == 12
        h = 10;
    end
    if flags(partId) == 0
        continue;
    end
    part = parts{partId};
    if abs(part.frame.center(1)) < 3e-2
        parts_non_pairs{length(parts_non_pairs)+1} = part;
    else
        tp = 0;
        for j = (partId+1):numParts
            part1 = parts{j};
            d1 = part1.frame.center(1) + part.frame.center(1);
            d2 = part1.frame.center(3) - part.frame.center(3);
            d3 = part1.frame.center(2) - part.frame.center(2);
            if abs(d1) < 3e-2 && abs(d2) < 3e-2 & abs(d3) < 3e-2
                flags(j) = 0;
                if j == 12
                    h = 10;
                end
                parts_pairs{length(parts_pairs)+1} = part;
                tp = 1;
                break;
            end
        end
        if tp == 0
            parts_non_pairs{length(parts_non_pairs)+1} = part;
        end
    end
end

%for i = 1:length(parts_pairs)
%  save2obj(parts_pairs{i}.mesh, sprintf('E:\\%d.obj', i));
%end

%for i = 1:length(parts_non_pairs)
%  save2obj(parts_non_pairs{i}.mesh, sprintf('E:\\%d.obj', i+length(parts_pairs)));
%end


function [frame] = extract_shape_part(Shape, mesh)
% Compute the frame of a shape part
    
numVertices = size(shape.vertexPoss, 2);
flags = zeros(1, numVertices);

flags(shape.faceVIds(:, faceIds)) = 1;
ids = find(flags > 0);
flags(ids) = 1:length(ids);

mesh.vertexPoss = shape.vertexPoss(:, ids);
mesh.faceVIds = flags(shape.faceVIds(:,faceIds));
[frame] = sp_comp_part_frame(mesh);

nv = size(mesh.vertexPoss, 2);
mesh.vertexPoss = frame.axis'*(mesh.vertexPoss - frame.center*single(ones(1,nv)));

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