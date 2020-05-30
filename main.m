clear;clc
p = genpath('devkit/');
addpath(p)

base_dir = 'D:\data\2011_10_03\2011_10_03_drive_0027_sync';
oxts = loadOxtsliteData(base_dir);

imu2velo = loadCalibrationRigid([base_dir, '\calib_imu_to_velo.txt']);
velo2cam = loadCalibrationRigid([base_dir, '\calib_velo_to_cam.txt']);
cam2imu = inv(velo2cam * imu2velo);

pos = getLocation(oxts);
m = imuMeasure(oxts);

mu_init = eye(4);
sigma_init = eye(15);
v_init = [oxts{1}(9);oxts{1}(10);oxts{1}(11)];
filter = LI_EKF(mu_init, sigma_init, v_init);

traj = [];
traj_gt = [];
for i = 2:length(pos)
    u = m(:, i);
    filter.prediction(u);
    p = pos(:, i);
    filter.correction(p);
    
    T = zeros(3, 4);
    T(1:3, 1:3) = filter.mu(1:3, 1:3);
    T(1:3, 4) = filter.mu(1:3, 5);
    traj(:, i) = reshape(T', [12, 1]);
    
    T = [T; [0, 0, 0, 1]];
    T = inv(cam2imu)*T;
    T = T(1:3, :);
    traj_gt(:, i) = reshape(T', [12, 1]);
end

dlmwrite('EKF_traj.txt', traj_gt', 'delimiter', ' ');
gt = readgt();

plot3(gt(1,:), gt(2, :), gt(3, :));
hold on
plot3(traj_gt(4, :), traj_gt(8, :), traj_gt(12, :));
xlabel('x')
ylabel('y')
axis equal
legend('Ground Truth', 'IEKF Trajectory')

function p = getLocation(oxts)
%     p = [0;0;0];
    rx = oxts{1}(4); % roll
    ry = oxts{1}(5); % pitch
    rz = oxts{1}(6); % heading 
    Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)]; % base => nav  (level oxts => rotated oxts)
    Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)]; % base => nav  (level oxts => rotated oxts)
    Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1]; % base => nav  (level oxts => rotated oxts)
    R  = Rz*Ry*Rx;
    for i = 1:length(oxts)
        [x, y, z] = geodetic2enu(oxts{i}(1),oxts{i}(2),oxts{i}(3),...
                                 oxts{1}(1),oxts{1}(2),oxts{1}(3),...
                                 wgs84Ellipsoid);
        
        p(:, i) = R'*[x;y;z];
    end
end

function m = imuMeasure(oxts)
    m = [];
    for i = 1:length(oxts)
        m(:, i) = [oxts{i}(18);oxts{i}(19);oxts{i}(20);
                   oxts{i}(12);oxts{i}(13);oxts{i}(14)];
    end
end