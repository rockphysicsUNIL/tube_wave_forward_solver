function [psi,grad_psi] = dispersion_relation_grad(kr,omega,L0,alpha_f,nu)
% Dispersion relation according to equation 4 of Bakku et al. 2013 and its gradient. 
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

% Calculation of the gradient
deta1_drealkr = -1.0.*(omega.^2./(alpha_f.^2-4.0./3.0.*1i.*omega.*nu)-(real_kr+1i.*imag_kr).^2).^(-1/2).*(real_kr+1i.*imag_kr);
deta1_dimagkr = -1.0.*(omega.^2./(alpha_f.^2-4.0./3.0.*1i.*omega.*nu)-(real_kr+1i.*imag_kr).^2).^(-1/2).*(real_kr+1i.*imag_kr).*1i;
deta2_drealkr = -1.0.*(1i.*omega./nu-(real_kr+1i.*imag_kr).^2).^(-1/2).*(real_kr+1i.*imag_kr);
deta2_dimagkr = -1.0.*(1i.*omega./nu-(real_kr+1i.*imag_kr).^2).^(-1/2).*(real_kr+1i.*imag_kr).*1i;
deta1eta2_drealkr = deta1_drealkr.*eta2 + eta1.*deta2_drealkr;
deta1eta2_dimagkr = deta1_dimagkr.*eta2 + eta1.*deta2_dimagkr;
df1_drealkr =     2.0.*(real_kr+1i.*imag_kr).*tan(eta2.*L0./2) +(real_kr+1i.*imag_kr).^2.*(cos(eta2.*L0./2)).^(-2).*L0./2.*deta2_drealkr;
df1_dimagkr = 1i.*2.0.*(real_kr+1i.*imag_kr).*tan(eta2.*L0./2) +(real_kr+1i.*imag_kr).^2.*(cos(eta2.*L0./2)).^(-2).*L0./2.*deta2_dimagkr;
df2_drealkr = deta1eta2_drealkr.*tan(eta1.*L0./2) +eta1.*eta2.*(cos(eta1.*L0./2)).^(-2).*L0./2.*deta1_drealkr;
df2_dimagkr = deta1eta2_dimagkr.*tan(eta1.*L0./2) +eta1.*eta2.*(cos(eta1.*L0./2)).^(-2).*L0./2.*deta1_dimagkr;
grad_psi(1,1) = 2.0.*real(f1+f2).*real(df1_drealkr+df2_drealkr) +2.0.*imag(f1+f2).*imag(df1_drealkr+df2_drealkr);
grad_psi(2,1) = 2.0.*real(f1+f2).*real(df1_dimagkr+df2_dimagkr) +2.0.*imag(f1+f2).*imag(df1_dimagkr+df2_dimagkr);
