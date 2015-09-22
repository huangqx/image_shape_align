function [Para] = set_default_parameters(Para_in)
Para = Para_in;

%Camera sampling
Para.numAzimuth = 36;
Para.numElevation = 5;
Para.elevationMin = 0.2000;
Para.elevationMax = 0.8000;
Para.lookAt = [0,0,0]';
Para.viewDistance = 3;
Para.scale = 0.4000;
Para.nHeight = 600;
Para.nWidth = 600;
Para.nHeight_inner = 734;
Para.nWidth_inner = 771;

% Hog computation
Para.gridHog = [6,6];
Para.numOrients = 16;

% Image-shape alignment
Para.knn_unary = 4;
Para.ffd_Res = 0.4;
Para.numIterations = 40;
Para.lambda_saliency = 0.1;
Para.minSigmaCorres = 0.001;
Para.lambda_1st = 0.5;
Para.lambda_smooth = 8;
Para.numIterations_non_rigid = 10;
Para.numIterations_rigid = 10;