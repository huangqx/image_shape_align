This software implements a simple image-shape alignment algorithm
% Load a dataset
load('data\car.mat');
% Each dataset consists of a few test images and a few shapes
% The shapes are consisttly oriented: the y-axis is the up-right orientation, and the underlying reflectional plane should be the yz-plane, the centers of the shapes are [0,0,0]
% Each image is given an estimated bounding box (does not need to be very accurate)

% Load the parameters
load('data\parameters.mat');

% Perform image-shape alignment for one image, it also select a cloest shape
[Camera_opt, Shape_opt] = demo_i2s_align(Images{1}, Shapes, Para);
