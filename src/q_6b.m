close all;
clear all;
clc;

%% PARAMETERS

% coefficients standing wave wavelength 1/2 lamprey
coeffs = [ -1.956393611785733e-06
            1.369475528250095e-04    
           -3.564099499976167e-03     
            4.125742217379403e-02    
           -1.767814629335283e-01
           -2.437207555102614e-01     
            3.163381532248638e+00    
           -4.028165285407771e+00]';

% linear drift propagating wave 2*pi phase difference at extremities
e = -0.092;

%% DATA
sw = chain_phase_oscil(10,-0.1,20);
pw = chain_phase_oscil(10,e,20);

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

