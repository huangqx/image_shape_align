function [Shape_out] = mm_color_shape(Shape_in, rgbd)
%
Shape_out = Shape_in;
field = 'crl_diffuse';
for i = 1:size(Shape_out.sym_pairs, 2)
    id1 = Shape_out.sym_pairs(1, i);
    id2 = Shape_out.sym_pairs(2, i);
    off = max(1, floor(rand(1,1)*64));
    rgb = rgbd(off, :);
    Shape_out.meshes{id1}.mat.clr_diffuse(1:3) = rgb;
    Shape_out.meshes{id2}.mat.clr_diffuse(1:3) = rgb;
end

for i = 1:length(Shape_out.nonsym_ids)
    id = Shape_out.nonsym_ids(i);
    off = max(1, floor(rand(1,1)*64));
    rgb = rgbd(off, :);
    Shape_out.meshes{id}.mat.clr_diffuse(1:3) = rgb;
end