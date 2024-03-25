clc
clear
N = 14; % (liczba segmentow - 1)
S = 5.1; % m^2
AR = 10; % Wydluzenie skrzydla
Czh=-0.202 %wsp silny nosnej usterzenia dla predkosci przelotowej 
lambda = 1; % Zbieznosc skrzydla
alpha_twist = 0.00000001; % Zwichrzenie geometryczne usterzenia(deg)
a_2d = 6.7/(1+6.7/(3.14*AR)); % dcz/dalpha (1/rad)
%a_h = (Czh/a_2d)*57.3; % kat natarcia usterzenia (deg)
a_h=-2.55
alpha_0 = 0.000000001; % kat natarcia dla zerowej sily nosnej (deg)
b = sqrt(AR*S); % rozpietosc usterzenia (m)
MAC = S/b; % srednia cieciwa aerodynamiczna usterzenia (m)
Croot = 2*S/(b*(1+lambda)); % root chord (m)
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = a_h+alpha_twist:-alpha_twist/(N-1):a_h; % kąt natarcia każdego semgmentu
z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % srednia cieciwa aerodynamiczna natarcia kazdego segmentu (m)
mu = c * a_2d / (4 * b);
LHS = mu .* (alpha-alpha_0)/57.3; % Lewa strona równania
% rozwiazywanie N rownan by znalezc wspolczynniki A_i
for i=1:N
    for j=1:N
    B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
    end
end
A=B\transpose(LHS);
for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j = 1 : N
        sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end
CL = 4*b*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9) CL(10) CL(11) CL(12) CL(13) CL(14)];
y_s=[b/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9) z(10) z(11) z(12) z(13) z(14)];

CL_tail = pi * AR * A(1)