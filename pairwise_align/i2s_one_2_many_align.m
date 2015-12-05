function [Camera_opt, shapeFFDs_opt] = i2s_one_2_many_align(inputImage,...
    Shapes,... % the input shapes
    Camera_init,... % The initial camera configuration
    Para)
% This function implements aligning one input image with multiple shapes
% We assume that the input shapes are consistently oriented in a world
% coordiante system
% It is expected that this is better than single pair image-2-shape
% alignments
% Input arguments:
%       inputImage.im:   the input image
%       inputImage.bbox: the boundingbox
%       Shapes:          the shape collection that is consistently oriented
%       Camera_init:     the initial camera configuration in the coordinate
%                        system associated with the world coordinate system
%                        of the Shapes
%       Para.numIterations: the number of alternating iterations
%       Para.lambda_1st:    the positional constraint used in non-rigid
%                            alignment
%       Para.lambda_2nd:    the smoothenss constraint used in non-rigid
%                            alignment
%       Para.lambda_cycle:  the panelty in front of the shape
%                            correspondences

function [shapeFFDs] = i2s_