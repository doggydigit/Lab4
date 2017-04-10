close all;
clear all;
clc;

%% PARAMETERS
% number of oscillators
Ns = [10 20];   % must be even

% variation of E
A  = 2;
Es = [-0.2 -0.1 -0.01];

% variation of A
E  = -0.1;
As = [0.5 2 10];

%% DATA

for i = 1:length(Ns)
    n = Ns(i);
    
    theta0 = 2*pi*(rand(n,1)-0.5);
    
    M = [];
    for j = 1:length(Es)
        e = Es(j);
        M = [M chain_phase_oscil(A,e,n)];
    end
    save(sprintf('results/6.a/N%i_A_const.mat',n),'n','A','Es','M');
    
    M = [];
    for j = 1:length(As)
        a = As(j);
        M = [M chain_phase_oscil(a,E,n)];
    end
    save(sprintf('results/6.a/N%i_E_const.mat',n),'n','As','E','M');
end

%% PLOTS
close all;
clc;

% Arnold tongue plot
figure(1);
hold on;
grid on;

xlabel('Drift [rad]');
ylabel('Coupling strength');
xlim([-0.25 0.25]);
ylim([0 15]);

xL = xlim;
yL = ylim;
area([xL(1) 0 xL(2)],[abs(Ns(1)^2*xL(1)/8) 0 abs(Ns(1)^2*xL(2)/8)],yL(2),'LineStyle','--','FaceColor',[1 0.8 0]);
area([xL(1) 0 xL(2)],[abs(Ns(2)^2*xL(1)/8) 0 abs(Ns(2)^2*xL(2)/8)],yL(2),'LineStyle','-','FaceColor',[1 0.4 0]);
alpha(0.35);

plot([E E E],[As],'bo');
plot([Es],[A A A],'r*');

legend(sprintf('N = %i',Ns(1)),sprintf('N = %i',Ns(2)),sprintf('e = %.2f / a = [%.2f ; %.2f ; %.2f]',E,As),sprintf('a = %.2f / e = [%.2f ; %.2f ; %.2f]',A,Es));

% PHASE DIFFERENCE STD

for i = 1:length(Ns)
    load(sprintf('results/6.a/N%i_A_const.mat',Ns(i)));
    figure(2+2*(i-1));
    hold on;
    grid on;
    plot((1:size(M,1))*0.01,M);
    legend(sprintf('e = %.2f',Es(1)),sprintf('e = %.2f',Es(2)),sprintf('e = %.2f',Es(3)));
    xlabel('Time');
    ylabel('STD(Phase difference)');
    %savefig(sprintf('results/6.a/N%i_A_const.fig',Ns(i)));
    
    load(sprintf('results/6.a/N%i_E_const.mat',Ns(i)));
    figure(3+2*(i-1));
    hold on;
    grid on;
    plot((1:size(M,1))*0.01,M);
    legend(sprintf('a = %.2f',As(1)),sprintf('a = %.2f',As(2)),sprintf('a = %.2f',As(3)));
    xlabel('Time');
    ylabel('STD(Phase difference)');
    %savefig(sprintf('results/6.a/N%i_E_const.fig',Ns(i)));
end