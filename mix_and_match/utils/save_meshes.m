function [] = save_meshes(Shape, foldername)
%
for i = 1:length(Shape.meshes)
    mesh = Shape.meshes{i};
    mesh.vertexPoss = Shape.vertexPoss(:, mesh.vertexIds);
    save2obj(mesh, [foldername, sprintf('%d.obj', i)]);
end