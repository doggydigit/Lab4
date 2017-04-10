%% PLOTS
close all;
clc;

figure(1)
hold on;
grid on;

plot(-pw{1},'-o')

ylabel('Stable phase difference');
xlabel('Oscilators');

yyaxis right;
plot(cumsum(-pw{1})/(2*pi));
ylabel({'Stable phase difference with oscillator 1';'[fraction of 2\pi]'});

figure(2)
hold on;
grid on;

plot(-sw{1},'-o')

ylabel('Stable phase difference');
xlabel('Oscilators');

yyaxis right;
plot(-cumsum(sw{1})/(2*pi));
ylabel({'Stable phase difference with oscillator 1';'[fraction of 2\pi]'});