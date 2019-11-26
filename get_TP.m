function Tp = get_TP(mu,Ks,theta1,rho,eta,freq)
% Use the equations provided by Schoenberg (1980) to calculate the angle-dependent P-wave transmission coefficient. 
% Written by Jürg Hunziker, University of Lausanne
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etaT = eta; % Fracture compliance in transverse direction
etaN = eta; % Fracture compliance in normal direction
rho1 = rho; % Density above fracture
rho2 = rho; % Density below fracture
alpha1 = sqrt((Ks+4/3.*mu)./rho1); % P-wave velocity above fracture
alpha2 = sqrt((Ks+4/3.*mu)./rho2); % P-wave velocity below fracture
beta1 = sqrt(mu./rho1); % S-wave velocity above fracture
beta2 = sqrt(mu./rho2); % S-wave velocity below fracture
omega = 2*pi*freq; % angular frequency

theta2 = asind(sind(theta1)/alpha1*alpha2);
phi1 = asind(sind(theta1)/alpha1*beta1);
phi2 = asind(sind(theta1)/alpha1*beta2);
gamma1 = 2*rho1*beta1*sind(phi1);
gamma2 = 2*rho2*beta2*sind(phi2);
p1 = rho1*alpha1-gamma1*sind(theta1);
p2 = rho2*alpha2-gamma2*sind(theta2);
q1 = rho1*beta1*(cosd(phi1))^2-0.5*gamma1*sind(phi1);
q2 = rho2*beta2*(cosd(phi2))^2-0.5*gamma2*sind(phi2);

A = [-p1, gamma1*cosd(phi1), p2, gamma2*cosd(phi2);...
     gamma1*cosd(theta1), q1, gamma2*cosd(theta2), -q2;...
     -sind(theta1), -cosd(phi1), sind(theta2)-i*omega*etaT*gamma2*cosd(theta2), -cosd(phi2)+i*omega*etaT*q2;...
     cosd(theta1), -sind(phi1), cosd(theta2)-i*omega*etaN*p2, sind(phi2)-i*omega*etaN*gamma2*cosd(phi2)]; 
Bp = [-A(1,1);A(2,1);-A(3,1);A(4,1)];

x = inv(A.'*A+1e-6)*A.'*Bp;
Tp = x(3);
