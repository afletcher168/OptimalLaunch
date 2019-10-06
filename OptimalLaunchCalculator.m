%%File: OptimalLaunchCalculator.m
%%Author: Andrew Fletcher
%%This prgram can be used for a few purposes pertaining to launch
%%characteristics. The probably most useful function is it's ability to
%%take a launch distance and height, and calculate the minimum launch speed
%%needed (by launching at the optimal angle)(also returns the launch angle)
%%Due to the unknown nature of the object (it's shape/mass/surface
%%characteristics and the characteristics of the air), air resistance is
%%neglected.

%d/dx = v0*cosd(2*theta) - sind(theta)*sqrt(v0^2*(sind(theta))^2 + 2*g*h0) + v0^2*sind(theta)*(cosd(theta))^2/sqrt(v0^2*(sind(theta))^2 + 2*g*h0);
%when d/dx = 0, theta = optimal theta, but that is a very hard equation to
%solve. To loop through all the theta values to find when d/dx is around 0
%would be redundant with what this program already does.
%According to a website answer, the root to d/dx = 0 is:
%theta = acos(sqrt(2*g*h0/(v0^2) +1)/sqrt(2*g*h0/(v0^2) + 2))
%This formula seems to work for finding the optimal theta.

clear;
clc;
close all;
format shortg;

MAX_RANGE = false; %set to true to find angle for max range, false for min speed to land at set distance

%%NOTE: MAKE SURE TO KEEP UNITS CONSISTENT THROUGHOUT PROGRAM
g = 9.81;  %gravitational acceleration constant on earth
h0 = -4;  %initial height
%theta = 45;    %launch angle
if MAX_RANGE
    v0 = 30;   %specified initial speed
else
    xDistance = 22.25;    %x distance of landing point
    v0 = 0.1:0.1:100; %checking every initial speed 0.1-100 with step = 0.1
end
theta = 0:0.1:90;    %checking every launch angle 0-90 with steps of 0.1

travelTime = zeros(length(v0),length(theta));
distance = zeros(length(v0),length(theta));

maxRange = zeros(1,length(v0));
optTheta = zeros(1,length(v0));
minSpeed = zeros(1,length(v0));
v0Counter = 0;
for v0 = 0.1:0.1:100
    v0Counter = v0Counter + 1;
    thetaCounter = 0;
    for theta = 0:0.1:90
        thetaCounter = thetaCounter + 1;
        if (v0^2 * (sind(theta))^2) < -2*g*h0
            travelTime(v0Counter,thetaCounter) = 0;
            distance(v0Counter,thetaCounter) = 0;
            continue
        end
        travelTime(v0Counter,thetaCounter) = (v0*sind(theta) + sqrt(v0^2 * (sind(theta))^2 + 2*g*h0))/g;
        distance(v0Counter,thetaCounter) = travelTime(v0Counter,thetaCounter) * v0 * cosd(theta);
        if distance(v0Counter,thetaCounter) > maxRange(v0Counter)
            maxRange(v0Counter) = distance(v0Counter, thetaCounter);
            optTheta(v0Counter) = theta;
            minSpeed(v0Counter) = v0;
        end
    end
end

optimalLaunch = [maxRange; optTheta; minSpeed];
theta = 0:0.1:90;

if MAX_RANGE
%%FOR A GIVEN LAUNCH SPEED
%Mathematically calculated optimal theta
formulaOptTheta = acosd(sqrt(2*g*h0/(v0^2) +1)/sqrt(2*g*h0/(v0^2) + 2));
sprintf("For a given speeed of %.2f, the maxmimum range occurs at %.1f degrees.(as calculated by formula)",...
    v0, formulaOptTheta)

%%Finding max range at a GIVEN LAUNCH SPEED and height at optimal launch angle
sprintf("For a given speeed of %.2f, the maxmimum range is %.2f at %.1f degrees.",...
    v0, maxRange, optTheta)

%Plotting range vs. theta for GIVEN LAUNCH SPEED
figure(4)
plot(theta, distance, 'r')
xlabel("range");ylabel("launch angle (degrees");legend("theta")
title("range vs. launch angle for given speed");

else
%%FINDING MINIMUM LAUNCH SPEED
v0 = 0.1:0.1:100;
%%Plotting optimal launch angle
figure(1)
index1 = find(abs(v0-20) < 0.001,1);  %%Finding the index of when launch speed = 20
index2 = find(abs(v0-40) < 0.001,1);
index3 = find(abs(v0-60) < 0.001,1);
index4 = find(abs(v0-80) < 0.001,1);
index5 = find(abs(v0-100) < 0.001,1);
plot(theta,distance(index1,:),'r', theta,distance(index2,:),'m', theta,distance(index3,:),'b', theta,distance(index4,:),'c', theta,distance(index5,:),'g');
xlabel('theta (degrees)');ylabel('distance');legend('v0=20','v0=40','v0=60','v0=80','v0=100');
title('distance vs theta for various v0');
 
figure(2)
plot(maxRange,optTheta,'r');
xlabel('max range');ylabel('optimal theta (degrees)');legend('optimal theta');
title('optimal angle vs range');

figure(3)
plot(maxRange,minSpeed,'b')
xlabel("range");ylabel("minimum launch speed");legend("launch speed");
title("launch speed vs range at optimal angle");

%%Finding min power at optimal launch angle for given distance
index = find(maxRange >= xDistance,1);
sprintf("For travelling %.2f units, from a height of %.2f, the minimum speed is %.1f at %.1f degrees.",...
        xDistance, h0, minSpeed(index), optTheta(index))
end





