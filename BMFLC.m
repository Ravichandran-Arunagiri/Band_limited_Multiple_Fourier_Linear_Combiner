t = 0 : 0.001 : 120;

Mpv = sin(t)+2*sin(2*t)+3*sin(3*t)+ cos(t)+2*cos(2*t)+3*cos(3*t); % voluntary response
Mpi =  0.1*sin(30*t)+0.1*sin(35*t)+ 0.1*cos(30*t)+ 0.1*cos(35*t); % involuntary component
Mp = Mpv+Mpi;
x = Mp;


%%
%define minimum and max frequencies
f_min = 0;
delta_f = 0.01;
f_max = 50;

%mu - tuning factor
mu = 0.000003;

harmonics = 1; %beta
L = harmonics*f_max+1;
N = L/delta_f;
y_ref = x;
v =zeros(N,1);
z = zeros(2,N);
w = zeros(2,N);

%create an array containing all freqencies
for i = 1:N 
    v(i) =f_min+delta_f*i/harmonics;
end


y = zeros(length(t),1);
E = zeros(length(t),1);

for i = 1: length(t)
    for j =1:N
        z(1,j)=sin(v(j)*t(i)); 
        z(2,j)=cos(v(j)*t(i));
        
    end
    for j =1:N
        E(i) = y_ref(i) - y(i);
        w(1,j) = w(1,j)+2*mu*z(1,j)*E(i);
        w(2,j) = w(2,j)+2*mu*z(2,j)*E(i);
    end
    for j =1:N
        y(i) = y(i) + w(1,j)*z(1,j) + w(2,j)*z(2,j); 
    end
    
end
theta = [w(1,:)';w(2,:)'];


[m n] = size(theta);
f_av = 0;
f_bv = 5;
f_ai = 25;
f_bi = 50;

gamma = zeros(8,1);
gamma(1) = harmonics*(f_av-f_min)/delta_f+1;
gamma(2) = harmonics*(f_bv-f_min)/delta_f+1;
gamma(3) = L/delta_f+gamma(1);
gamma(4) = L/delta_f+gamma(2);
gamma(5) = harmonics*(f_ai-f_min)/delta_f+1;
gamma(6) = harmonics*(f_bi-f_min)/delta_f+1;
gamma(7) = L/delta_f+gamma(5);
gamma(8) = L/delta_f+gamma(6);

theta_v = [theta(gamma(1):gamma(2));theta(gamma(3):gamma(4))];
theta_i = [theta(gamma(5):gamma(6));theta(gamma(7):gamma(8))];
y_v = zeros(length(t),1);
y_i = zeros(length(t),1);
for i = 1: length(t)
     for j =1:N
        z(1,j)=sin(v(j)*t(i));
        z(2,j)=cos(v(j)*t(i));
    end
    phi = [z(1,:)';z(1,:)'];
    phi_v = [phi(gamma(1):gamma(2));phi(gamma(3):gamma(4))];
    phi_i = [phi(gamma(5):gamma(6));phi(gamma(7):gamma(8))];
    y_v(i) = theta_v'*phi_v;
    y_i(i) = theta_i'*phi_i;
end
close all
figure()
plot(t,y_ref,'r')
hold on
plot(t,y,'b')
legend('actual','estimated')
title('Mp')
hold off

figure()
plot(t,Mpv,'r')
hold on
plot(t,y_v,'b')
legend('actual','estimated')
title('Mpv')
hold off







