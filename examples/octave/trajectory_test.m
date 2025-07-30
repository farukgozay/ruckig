clc; clear all; close all;
clear sCurveRealtime1; 

%Parameters
CurrentPos  = 0;     
TargetPos1  = 20;    
TargetPos2  = -10;     
TargetPos3  = 45;      
TargetPos4  = 0;      
v_max       = 30;       
a_max       = 25;       
d_max       = 12;       
jerk_max    = 1000;      
Ts          = 0.01;     


time = 0:Ts:50;
N = length(time);


position     = zeros(N,1);
velocity     = zeros(N,1);
acceleration = zeros(N,1);
deceleration = zeros(N,1);  
jerk         = zeros(N,1);



TargetPos     = TargetPos1;
TargetPosPrev = TargetPos;
prev_vel=0;

for i = 1:N
   
    if time(i) >= 6
        TargetPos = TargetPos4;
    elseif time(i) >= 4
        TargetPos = TargetPos3;
    elseif time(i) >= 2
        TargetPos = TargetPos2;
    end
    
    if TargetPos ~= TargetPosPrev

        resetFlag = true;
        TargetPosPrev = TargetPos;
        CurrentPos = position(i-1);
        v_init = velocity(i-1);
        a_init = acceleration(i-1);
    else
        resetFlag = false;
    end




    [pos, vel, acc, dec, jerk_i] = sCurveRealtime1(TargetPos, CurrentPos, Ts, v_max, a_max, d_max, jerk_max);

    if abs(vel - prev_vel) > v_max
    vel =prev_vel + sign(vel- prev_vel) * v_max* Ts;
    end
    prev_vel=vel;



    position(i)     = pos;
    velocity(i)     = vel;
    acceleration(i) = acc;
    deceleration(i) = dec; 
    jerk(i)         = jerk_i;


    CurrentPos = pos;
end
figure('Name','Adaptif S-Curve Trajectory','NumberTitle','off');

subplot(4,1,1);
plot(time, position, 'b', 'LineWidth', 1.5); hold on;
grid on; ylabel('Position');
title('S-Curve Position');

subplot(4,1,2);
plot(time, velocity, 'r', 'LineWidth', 1.5);
ylim([-v_max*1.2, v_max*1.2]);   % ← burası eklendi
grid on; ylabel('Velocity');

subplot(4,1,3);
plot(time, acceleration, 'g', 'LineWidth', 1.5);
grid on; ylabel('Acceleration');


subplot(4,1,4); hold on;
plot(time, jerk, 'm', 'LineWidth', 1.5);
grid on; ylabel('Jerk'); xlabel('Time (s)');
ylim([-jerk_max*1.2, jerk_max*1.2]);
