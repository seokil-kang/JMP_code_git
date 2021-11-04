clear;
close all;
clc;
prif = [4, 1, 0.5];
prif = [prif;4, 1.2, 0.5];
N = 1e5;

x = prior_drawing(prif,N);
alph = 0;
quantile(x,[alph .5 1-alph],2)
%%
figure
lw = 2;
hold on
histogram(x(1,:),'normalization','pdf','displaystyle','stairs','linewidth',lw);
histogram(x(end,:),'normalization','pdf','displaystyle','stairs','linewidth',lw);
hold off