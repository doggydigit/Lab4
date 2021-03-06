function [ synchro ] = ex7bsimul(a, f)
% Two coupled neural oscillators made of leaky-integrator neurons
global tau;
global D;
global b;
global w;
global a_sens;
global omega_mech;

% Parameters of the neural network
% MODIFY THE PARAMETER VALUES

intfreq = 0.04;

tau   = [intfreq, intfreq, intfreq*f, intfreq*f]';  % CHANGE THIS TO CHANGE THE INTRINSIC FREQUENCIES
D     = 1;

b     = [-2.75, -1.75, -2.75, -1.75]';     % Values for a limit cycle
w     = zeros(4,4);
w(1:2,1:2)     = [4.5, -1; 1, 4.5]';  % First oscillator
w(3:4,3:4)     = [4.5, -1; 1, 4.5]';  % Second oscillator

% Coupling between the two neural oscillators:
% ADD COUPLING WEIGHTS HERE BETWEEN NEURONS 1 AND 3 (OR OTHER NEURONS OF YOUR CHOICE).
w(1,3) = a;
w(3,1) = a;

% Forcing term from a periodic mechanical input
% SET VALUES FOR THE PERIODIC FORCING TERM, see question 7.b
a_sens = 0;
omega_mech = 2*pi*0.8;

% Initial conditions
%  SET THE INITIAL CONDITIONS
y_0 = [0 1 2 3]';   % Values of the membrane potentials of the four neurons

dt = 0.001; % Force ode45 to return values at small steps
% RungeKutta Integration
[T,Y] = ode45(@(t,y) LI_network_ode(t,y),[0:dt:20],y_0);

% Recompute the firing rates for the figures:
X(:,1) = 1./(1+exp(-D*(Y(:,1)+b(1))));
X(:,2) = 1./(1+exp(-D*(Y(:,2)+b(2))));
X(:,3) = 1./(1+exp(-D*(Y(:,3)+b(3))));
X(:,4) = 1./(1+exp(-D*(Y(:,4)+b(4))));

% Compute the phases using the Hilbert transform:
z1 = hilbert(Y(:,1)-mean(Y(:,1)));
phase1 = angle(z1);
unwrapped_phase1 = unwrap(angle(z1));
% Estimate the instantaneous and average frequencies (ignore beginning and end of signal)
ignored_steps = 2000;
freq1 = (unwrapped_phase1(end-ignored_steps)-unwrapped_phase1(ignored_steps)) /  (T(end-ignored_steps)-T(ignored_steps));
freq1=freq1/(2*pi);

z3 = hilbert(Y(:,3)-mean(Y(:,3)));
phase3 = angle(z3);
unwrapped_phase3 = unwrap(angle(z3));
% Estimate the instantaneous and average frequencies (ignore beginning and end of signal)
freq3 = (unwrapped_phase3(end-ignored_steps)-unwrapped_phase3(ignored_steps)) /  (T(end-ignored_steps)-T(ignored_steps));
freq3=freq3/(2*pi);

if abs(freq3-freq1)<0.02
    synchro = 1;
else
    synchro = 0;
end



% Make a figure of the time evolutions of m and x
%close all;
% figure(1)
% subplot(2,1,1);
% axis([0 T(end) 0 5])
% set(gca,'FontSize',20)
% hold on
% plot(T,Y(:,1:2),'LineWidth',2)
% legend('Neuron1','Neuron2')
% xlabel('time'); ylabel('m');
% subplot(2,1,2);
% axis([0 T(end) 0 5])
% set(gca,'FontSize',20)
% hold on
% plot(T,Y(:,3:4),'LineWidth',2)
% legend('Neuron3','Neuron4')
% xlabel('time'); ylabel('m');



end


function yd = LI_network_ode(t,y)
% Derivative function of a network of leaky integrator neurons,
% y is the vector of membrane potentials (variable m in lecture equations)

global tau;
global D;
global b;
global w;
global a_sens;
global omega_mech;

% update the firing rates:
x = 1./(1+exp(-D*(y+b)));

% compute the dentritic sums for both neurons
dend_sum =  w *x;

% Apply a periodic forcing term to one of the neural oscillators
% (mechanical input). Here it is applied to Neuron 1. FEEL FREE TO CHANGE
dend_sum(1) = dend_sum(1) + a_sens* sin(omega_mech*t);

% compute the membrane potential derivative:
yd = (dend_sum-y)./tau;

end


