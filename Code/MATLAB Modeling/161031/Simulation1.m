%%%%%%
%Simulation
clear all
clc
global puck mallet table fig score goalsize comp e;
global width height input output mutability crossover population genome;

puck.d = 1.5;
puck.m = .1;
mallet.m = 5;
mallet.d = 2;
mallet.v = [0,0];
comp.m = mallet.m;
comp.d = mallet.d;
comp.v = [0,0];
puck.v = [0,15];
score = [0,0];
e = 0.9;

goalsize = 4;

table.x = 12;
table.y = 24;

fig = figure('units','normalized','position',[0 0 .3 .8]);
title(sprintf('Human   %d - %d   Computer',score))
axis equal
hold on;
table.obj = rectangle('Position',[0,0,table.x,table.y]);
rectangle('Position',[table.x/2-goalsize/2,-1,goalsize,2],'FaceColor','w','EdgeColor','none');
rectangle('Position',[table.x/2-goalsize/2,table.y-1,goalsize,2],'FaceColor','w','EdgeColor','none');
puck.obj = rectangle('Curvature',[1,1],'Position',[3.9,10,puck.d,puck.d],'FaceColor','r');
mallet.obj = rectangle('Curvature',[1,1],'Position',[0,0,mallet.d,mallet.d],'FaceColor','b');
comp.obj = rectangle('Curvature',[1,1],'Position',[table.x/2 - mallet.d/2,table.y - comp.d - 3,comp.d,comp.d],'FaceColor','g');

%AI VARIABLES
width = 20;
height = 2;
input = 6;
output = 2;
mutability = 0.001;
crossover = 0.75;
population = 40;

play();


