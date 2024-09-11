% read fom and rom data
fx=dlmread('fom.dragx.dat');
fy=dlmread('fom.dragy.dat');
rx=dlmread('rom.dragx.dat');
ry=dlmread('rom.dragy.dat');

% plot fom result:

tstart=401;
tend=500;

figure;
plot(fx(:,2),fx(:,3)); xlim([tstart tend]);
title('Drag result');
xlabel('Time');
ylabel('Drag');

% plot rom result:
%figure;
hold on
rx(:,1)=rx(:,1)+tstart; % restart from t=400
%plot(rx(:,1),rx(:,4)); xlim([400 500]);
plot(rx(:,1),rx(:,4),'--');
legend('FOM','ROM');
%title('ROM drag result');
xlabel('Time');
ylabel('Drag');

