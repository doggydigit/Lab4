function chain_phase_oscil(a,e,Noscils)
 % a is the coupling strength. e is the difference of intrinsic frequencies
 % from one oscillator to the next (typically negative, to get a head to
 % tail traveling wave). Noscils is the number of oscillators.
 close all
if nargin<1
    a = 4; % Coupling strength. FEEL FREE TO CHANGE
end
if nargin<2
    e=-pi/6;  % Difference of intrinsic frequencies from one oscillator to the next. FEEL FREE TO CHANGE
end
if nargin<3
    Noscils=10;  % Number of oscillators. FEEL FREE TO CHANGE
end
close all

disp(sprintf('\n Noscils: %d,   Difference of freq: e=%.2f  Coupling strength a=%0.2f\n',Noscils,e,a))

% Set the frequency values to the oscillators. CHANGE FOR QUESTIONS 6.b
% Oscillator 1 is the head. Oscillator Noscils is the tail
omega = 2*pi*ones(Noscils, 1);
for i=2:Noscils
    omega(i)=omega(i-1)+e;
    %omega(i)=omega(1)+exp((i-1/Noscils)^5)*e*Noscils;
    %omega(i)=omega(i-1)*0.96;
end
dif_omega = diff(omega);


% USE HERE WHAT YOU KNOW FROM THE ANALYSIS PRESENTED IN THE LECTURE TO
% PREDICT WHETHER THE SYSTEM WILL CONVERGE OR NOT.

% COMPUTE THE COUPLING MATRIX
d = -2*ones(Noscils-1,1);
d1= ones(Noscils-2,1);
A = a*(diag(d,0) + diag(d1,1) + diag(d1,-1));

% Compute the fixed points
% COMPUTE THE VECTOR S
S = -(A\dif_omega);

% DETERMINE THE SYNCHRONIZATION CRITERION: (example
% synchronization_achieved=a>12 )
synchronization_achieved = a>=abs(e)*Noscils*Noscils/8;

if synchronization_achieved
    disp(sprintf('Synchronization achieved:\n'));
    % COMPUTE THE FIXED POINTS FOR THE PHASE DIFFERENCES:
    fixed_points = asin(S);
        
    figure(1)
    plot(fixed_points,'o-','LineWidth',2)
    set(gca,'FontSize',20)
    set(gca,'YDir','Reverse')
    aa=axis;
    axis([0 Noscils aa(3) 0])
    xlabel('Oscillator')
    ylabel('Stable phase difference')
    hold off
    print -dpng chain_phase_oscil_phases.png
end

% Running the dynamical system

dt = 0.01;
t_end = 20;

% evolution of the different oscillators
dtheta = zeros(Noscils,1);
theta  = 2*pi*(rand(Noscils,1)-0.5); % Random initial conditions

counter = 0;
% Euler integration (NOTE: FEEL FREE TO IMPLEMENT AS A RUNGE-KUTA
% INTEGRATION)
for t=0:dt:t_end
    counter = counter+1;

        % apply a mechanical input to the tail segment
    % UPDATE HERE THE TERM TO APPLY A PERIODIC FORCING TERM, CF QUESTION
    % 6.C
    a_sens = 0;
    omega_mech = 2*pi*2;
    
    % Evolution of the different oscillators
    % COMPUTE THE DERIVATIVES OF THETA FOR THE DIFFERENT OSCILLATORS    
    dtheta = omega + [a*sin(theta(2:end)-theta(1:end-1)); a_sens*sin(omega_mech*t-theta(end))] + [0; a*sin(theta(1:end-1)-theta(2:end))]; 
   
    % Euler integration
    theta = theta + dtheta*dt;
    
    %keep logs
    T(counter) = t;
    THETA(counter,:) = theta;
    DTHETA(counter,:) = dtheta;
    X(counter,:)  = cos(theta);
end

% Estimate the resulting frequencies of each oscillators by taking the
% average rate of each oscillators for the last 20 time steps
for i=1:Noscils
    resulting_freq(i) = mean(DTHETA([end-20 end], i));
end

% Make figures and plots

figure(2)
subplot(211)
hold on
set(gca,'FontSize',20)
for (i=1:Noscils)
    X_shifted(:,i) = X(:,i) + 1 + Noscils-i;
end
plot(T,X_shifted,'LineWidth',2)
axis([0 t_end -0.2 Noscils+1.2])
set(gca,'ytick',[])
xlabel('Time')
ylabel('Oscillations')

subplot(212)
hold on
set(gca,'FontSize',20)
plot(T,mod(diff(THETA')',2*pi),'LineWidth',2)
if abs(S)<=1.0
    plot(T(end),fixed_points,'o')
end
axis([0 t_end 0 2*pi])
xlabel('Time')
ylabel('Phase differences')
print -dpng chain_phase_oscil.png

figure(3)
hold   on
set(gca,'FontSize',20)

if   a_sens == 0.0
    plot([1: Noscils],omega,'o-',[1: Noscils],resulting_freq,'o-','LineWidth',2)
    legend('Intrinsic', 'Resulting','location','best')
else
    plot([1: Noscils],omega,'o-',[1: Noscils],resulting_freq,'o-',Noscils,omega_mech,'o','LineWidth',2)
    legend('Intrinsic', 'Resulting','Mechanical','location','best')
end

axis([1 Noscils 0 max(omega)+10])
xlabel('Oscillator')
ylabel('Frequencies')
print -dpng chain_phase_oscil_freq.png
end
