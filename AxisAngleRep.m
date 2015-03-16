function [theta w] = AxisAngleRep(R)
%
% Convert a rotation matrix into axis-angle representation
%
% math from: http://en.wikipedia.org/wiki/Axis-angle_representation
%


theta = acos(.5*(trace(R)-1));

w = .5 * [ R(3,2) - R(2,3); ...
    R(1,3) - R(3,1);...
    R(2,1) - R(1,2)]./sin(theta);
