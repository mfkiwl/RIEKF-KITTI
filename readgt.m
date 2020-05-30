function gt = readgt()

addpath("C:\Users\Owner\Desktop\dataset\poses");
filename = '00.txt';
file = fopen(filename, 'r');

tline = fgetl(file);
traj = {};
while ischar(tline)
    arr = strsplit(tline,' ');
    pose = eye(4);
    for i = 1:12
        r = idivide(int32(i), 4, 'ceil');
        c = i-(r-1)*4;
        pose(r, c) = str2double(arr{i});
    end
    tline = fgetl(file);
    traj{length(traj)+1} = pose;
end
fclose(file);

gt = [];
for i = 1:length(traj)
    p = traj{i}(1:3, 4);
    gt = [gt, p];
end

end


