function [Camera_nor] = i2s_normalize_camera(object_bbox,...
    Shape, Camera_in)
%
% Modify the camera parameters so that the rendered objects align
% with the given object bounding boxes
render_image = i2s_render_shape(Shape, Camera_in);
render_bbox = i2s_object_bbox(render_image);
dims_o = [object_bbox(3) - object_bbox(1), object_bbox(4) - object_bbox(2)];
dims_r = [render_bbox(3) - render_bbox(1), render_bbox(4) - render_bbox(2)];
ratio = norm(dims_o)/norm(dims_r);
Camera_nor = Camera_in;
Camera_nor.scale = Camera_in.scale/ratio;

render_image = i2s_render_shape(Shape, Camera_nor);
render_bbox = i2s_object_bbox(render_image);
render_center = [(render_bbox(1)+render_bbox(3))/2,...
    (render_bbox(2)+render_bbox(4))/2];

object_center = [(object_bbox(1)+object_bbox(3))/2,...
    (object_bbox(2)+object_bbox(4))/2];

d = object_center - render_center;
[m, n, k] = size(render_image);
d = d*2*Camera_nor.scale/min(m,n);

axis_y = Camera_nor.upVec;
axis_z = Camera_nor.lookAt - Camera_nor.origin;
axis_z = axis_z/norm(axis_z);
axis_x = cross(axis_y, axis_z);
Camera_nor.lookAt = Camera_nor.lookAt + axis_x*d(1) + axis_y*d(2);
Camera_nor.origin = Camera_nor.origin + axis_x*d(1) + axis_y*d(2);







function [bbox] = i2s_object_bbox(render_image)

[m,n,k] = size(render_image);
if k > 1
    [Gmag, Gdir] = imgradient(rgb2gray(render_image), 'roberts');
else
    [Gmag, Gdir] = imgradient(render_image, 'roberts');
end
Gmag = abs(Gmag);
val = max(max(Gmag));
cols = max(Gmag);
cols = find(cols > val/100);
rows = max(Gmag');
rows = find(rows > val/100);
bbox = [min(cols), min(rows), max(cols), max(rows)];

