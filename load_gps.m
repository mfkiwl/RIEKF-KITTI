clear;clc
p = genpath('devkit/');
addpath(p)

base_dir = 'D:\data\2011_09_30_drive_0016_sync';
oxts = loadOxtsliteData(base_dir);
pose = convertOxtsToPose(oxts);

% plot every 10'th pose
figure; hold on; axis equal;
% for i=1:10:length(pose)
for i=1:10:20
  trans = pose{i}(1:3, 4);
  r = pose{i}(1:3, 1:3);
  rot = compact(quaternion(r, 'rotmat', 'frame'));
  plotTransforms(trans', rot, 'FrameSize', 10);
end
xlabel('x');
ylabel('y');
zlabel('z');

filter = RI_EKF();