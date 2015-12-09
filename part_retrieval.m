function [Shape_init, partCorres] = part_retrieval(Image, Shapes,...
    partCorresStructs, Para)
% This function retrieves a subset of parts of the input shapes to
% reconstruct the input image object. This is based on solving a
% generalized set-cover problem


