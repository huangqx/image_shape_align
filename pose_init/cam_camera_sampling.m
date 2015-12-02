function [Cameras] = cam_camera_sampling(Para)
% Para.numAzimuth = 60, number of samples for azimuth
% Para.numElevation = 6, number of samples for elevation
% Para.elevationMin = 0.2, the minimum camera elevation
% Para.elevationMax = 0.8, the maximum camera elevation
% Para.lookAt: the default lookAt position, which is usually
%              set at [0, 0.0, 0]'
% Para.viewDistance: the default viewing distnace, which is set
%                    at 3.2
% Para.scale: the default scaling parameter in matlab rendering, 
%             which is set to 0.4 (so that the model is included 
%             in the rendered image, yet it occupies the majority
%             of the image
% The following are default parameters of the camera model
% Para.nHeight = 600;
% Para.nWidth = 600;
% Para.nHeight_inner = 734;
% Para.nWidth_inner = 771;

for i = 1:Para.numElevation
    t = (i-1)/(Para.numElevation-1);
    phi = (1-t)*Para.elevationMin + t*Para.elevationMax;
    for j = 1:Para.numAzimuth
        cameraId = (i-1)*Para.numAzimuth + j;
        theta = 2*pi*(j-0.5)/Para.numAzimuth;
        Camera.lookAt = Para.lookAt;
        d = [sin(theta)*cos(phi), sin(phi), cos(theta)*cos(phi)]';
        Camera.origin = Camera.lookAt + d*Para.viewDistance;
        Camera.upVec = [0, 1, 0]';
        Camera.upVec = Camera.upVec - d*(Camera.upVec'*d);
        Camera.upVec = Camera.upVec/norm(Camera.upVec);
        Camera.scale = Para.scale;
        Camera.nHeight = Para.nHeight;
        Camera.nWidth = Para.nWidth;
        Camera.nHeight_inner = Para.nHeight_inner;
        Camera.nWidth_inner = Para.nWidth_inner;
        Cameras{cameraId} = Camera;
    end
end
