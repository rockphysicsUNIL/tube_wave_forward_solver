function kr_out = find_kr_CG(L0,Kf,rhof,nu,omega)
% Conjugate gradient algorithm to find the characteristic wavenumber kr. 
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
%
% Based on the paper: Shewchuk, J. R., 1994, An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
% https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf

% General Conjugate Gradient Inversion Parameters
whichstep = 2; % Determine the step length by (1) fitting of a parabola or 
               % (2) by a line search.
whichdir = 1; % Use (1) the Fletcher-Reeves formula to determine the search 
              % direction or (2) the Polak-Ribière formula. 
m = [0.05;0.05]; % starting guess for CG scheme [m^-1]
restartgrad = 3; % Amount of iterations after which the gradient is restarted.
ngradrestart = 4; % How many times an extraordinary restart of the gradient
                  % is allowed before the program stops.
maxit = 20; % maximum amount of iterations
eps = 10^(-6); % tolerance to stop search, i.e., the distance a new 
               % solution is away from the current solution: |alpha*d|
                 
% Line Search Parameters to determine the steplength in the CG scheme
maxita = 20; % maximum amount of iterations to find the step length alpha
epsa = 10^(-5); % tolerance to stop search for alpha
maxnogood = 5; % maximum amount of steps in a row that alpha has not improved

alpha_f = sqrt(Kf/rhof); % fluid velocity in the fracture [m/s]

% Using the conjugate gradient method
mvec = m;
ifun = 0; % This counts how many times the forward step has to be evaluated
[funm,grad] = dispersion_relation_grad(m,omega,L0,alpha_f,nu);
ifun = ifun + 1;
d = -grad;
r = -grad;

switcher = 0;
newgrad = 0;
in = 1;
lengthstep = 1;
stopcounter = 0;
while stopcounter<ngradrestart && in<=maxit && sum(grad)~=0
    if whichstep == 1 % Determining the steplength by fitting a parabola
        normfac = sqrt(sum(d.^2));
        kernel = [0 (1/normfac)^2 (2/normfac)^2;0 1/normfac 2/normfac;1 1 1];
        [funtemp1,gradtemp1] = dispersion_relation_grad([m(1)+1/normfac*d(1);m(2)+1/normfac*d(2)],omega,L0,alpha_f,nu);
        ifun = ifun + 1;
        [funtemp2,gradtemp2] = dispersion_relation_grad([m(1)+2/normfac*d(1);m(2)+2/normfac*d(2)],omega,L0,alpha_f,nu);
        ifun = ifun + 1;
        [detkernel,invkernel] = inv3x3(kernel);
        coef = [funm funtemp1 funtemp2]*invkernel;
        alpha = -coef(2)/(2*coef(1));
        if alpha<0
            % A negative alpha means that we are walking 180 degrees in the
            % wrong direction. 
            alpha = 1/normfac;
        end
        % fprintf('Parabola fitting: alpha = %e\n\n',alpha);
    elseif whichstep == 2 % Determining the steplength by a line search that uses only the function and not the gradient
        funtempprev = funm;
        alphafac = 1/(sqrt(sum(d.^2)));
        dalpha = 1;
        alpha = dalpha;
        alphaprev = alpha;
        inc = 0;
        dostop = 0;
        firsttry = 1;
        nogood = 0;
        ia = 0;
        while dostop == 0
            [funtemp,gradtemp] = dispersion_relation_grad([m(1)+alpha*alphafac*d(1);m(2)+alpha*alphafac*d(2)],omega,L0,alpha_f,nu);
            ifun = ifun + 1;
            % fprintf('Iteration %d: alpha = %e; J = %e; J0 = %e\n',ia,alpha,funtemp,funm);
            if funtemp<funtempprev
                % The new alpha leads to a smaller value of the function.
                % Alpha is slightly increased to see if that further
                % reduces the function value. 
                funtempprev = funtemp;
                alphaprev = alpha;
                while dalpha>alpha
                    % In case alpha has been previously reduced
                    % significantly, by entering firsttry==1 below, dalpha
                    % needs to be adjusted accordingly.
                    dalpha = dalpha/2;
                end
                alpha=alpha+dalpha;
                firsttry=0; % Alpha has been changed at least once this iteration.
                nogood=0; % The last change of alpha has reduced the functional. 
            else
                % The new alpha leads to a higher value in the function.
                if firsttry==1
                    % If firsttry==1 we have not updated alpha successfully
                    % before. Our first guess of alpha must have been to
                    % high. Alpha is therefore reduced significantly.
                    alphaprev=alphaprev/10;
                    alpha=alphaprev;
                else
                    % We have updated alpha successfully before. So, the
                    % first guess was good. We reduce dalpha to make a
                    % smaller new guess.
                    alpha = alphaprev;
                    dalpha = dalpha/2;
                    alpha=alpha+dalpha;
                    inc = 1; % The function has increased
                    nogood=nogood+1;
                end
            end
            if (inc==1 && ia>maxita) || (firsttry==1 && ia>maxita) || dalpha<epsa || nogood>maxnogood
                % We stop finding a better alpha if one of the following conditions is fulfilled:
                %
                % (inc==1 && ia>maxita) Alpha has increased, but we had a
                % good alpha before, because otherwise we would have ended
                % up in the firsttry==1 case and would not have set inc=1.
                % We have tried to find a new alpha for more then the
                % predefined tryouts in that case (maxita). 
                %
                % (firsttry==1 && ia>maxita)
                % Our first guess of alpha was not good, but we have
                % reduced that first guess more than the predefined tryouts
                % specified in maxita allow. 
                %
                % dalpha<epsa
                % The change in alpha, dalpha, is smaller than the predefined
                % value. 
                %
                % nogood>maxnogood
                % The variable nogood counts how many times in a row we
                % have tried a new alpha but only increased the function.
                % Nogood is reset to 0 if the function is reduced
                % intermittently. 
                dostop = 1;
                alpha = alphaprev; % The previous alpha is reinstated.
            end
            ia = ia + 1; % Counts the amount of new tryouts of alpha
        end
        % fprintf('------------------------------------------\n');
        % fprintf('Final alpha = %e\n\n',alpha);
        alpha = alpha*alphafac;
        if switcher == 1
            % Switcher is 1 if the previous iteration determined to switch
            % from parabola fitting to line-search. 
            whichstep = 1; % The line search is done and we switch back to parabola fitting
            switcher = 0; % reinitialize switcher
        end
    end
    
    mold = m;
    m = m+alpha*d;
    mvec(:,size(mvec,2)+1) = m;
    lengthstep = sqrt(sum((alpha*d).^2));
    funmold = funm;
    gradold = grad;
    [funm,grad] = dispersion_relation_grad([m(1);m(2)],omega,L0,alpha_f,nu);
    ifun = ifun + 1;
    if funm>funmold
        % fprintf('Function value for next solution has increased.\n')
        m = mold; % The previous parameters are reinstated
        funm = funmold;
        grad = gradold;
        if whichstep == 1 
            % If we were looking for alpha using parabola fitting, we next
            % time search for an alpha using the line search.
            % fprintf('Find a new alpha with a line search.\n')
            whichstep = 2; % do the line search next time
            switcher = 1; % Through this we switch back to parabola fitting after the line search
        end
    end
    % fprintf('Solution at Iteration %d: (%e,%e)\n',in,m(1),m(2));
    rold = r;
    rarray(in) = sqrt(sum(r.^2));
    r = -grad;
    if (lengthstep<eps && firsttry==1)
        % Because we were only in the fristtry==1 section, we have reduced
        % alpha to a very small amount. Therefore, lengthstep is smaller
        % than the predefined value. But that does not mean that we have
        % converged. It could also mean, that the search direction was not
        % optimal. Therefore, by setting newgrad=1, the search direction is
        % changed to the negative gradient for the next iteration and
        % lengthstep is boosted such that we don't stop iterating. 
        lengthstep=10*eps;
        newgrad = 1;
    end
    if lengthstep<eps
        % As above, lengthstep is very small, but not because we were in
        % the firsttry section. We could have been converged or the
        % direction is not good. The search direction is also changed to
        % the negative gradient but the counter to stop the whole process
        % is increased by 1 as well. 
        stopcounter = stopcounter + 1;
        newgrad = 1;
    end
    if mod(in,restartgrad)==0 || newgrad==1
        % Here the search direction is set to the negative gradient. This
        % happens if newgrad==1 but also at every ith iteration, which is 
        % specified by the variable restartgrad. 
        d = -grad;
        if newgrad == 1
            % The variable newgrad needs to be reinitialized, otherwise we
            % would end up here next iteration as well. 
            newgrad = 0;
        end
    else
        % In most iterations, i.e., when nothing special happened, we
        % compute the search direction with one of the two formulas below.
        if whichdir == 1
            beta = r.'*r/(rold.'*rold); % Fletcher-Reeves
        elseif whichdir == 2
            beta = r.'*(r-rold)/(rold.'*rold); % Polak-Ribière
        else
            error('Choose 1 or 2 for the parameter whichdir.');
        end
        d = r+beta*d;
    end
    
    in = in + 1; % Counts the amount of iterations
end

kr_out = complex(m(1),m(2));
