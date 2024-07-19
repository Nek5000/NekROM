% read fom and rom data
fx=dlmread('fom.dragx.dat');
fy=dlmread('fom.dragy.dat');
rx=dlmread('rom.dragx.dat');
ry=dlmread('rom.dragy.dat');

% plot fom result:
figure;
plot(fx(:,2),fx(:,3)); xlim([400 500]);
title('FOM drag result');
xlabel('Time');
ylabel('Drag');

% plot rom result:
figure;
plot(rx(:,1),rx(:,4)); xlim([400 500]);
title('ROM drag result');
xlabel('Time');
ylabel('Drag');

