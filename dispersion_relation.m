function psi = dispersion_relation(kr,omega,L0,alpha_f,nu)
% Dispersion relation according to equation 4 of Bakku et al. 2013.
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

real_kr = kr(1);
imag_kr = kr(2);

% Calculation of the objective function
eta1 = sqrt(omega.^2./(alpha_f.^2-4./3.*1i.*omega.*nu)-(real_kr+1i.*imag_kr).^2);
eta2 = sqrt(1i.*omega./nu-(real_kr+1i.*imag_kr).^2);
f1 = (real_kr+1i.*imag_kr).^2.*tan(eta2.*L0./2);
f2 = eta1.*eta2.*tan(eta1.*L0./2);
temp1 = (real(f1 + f2)).^2;
temp2 = (imag(f1 + f2)).^2;
psi = temp1+temp2;
