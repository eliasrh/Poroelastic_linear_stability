%% Example code for solving for the critical wavenumber kcrit 

% This code implements and constructs the charateristic equation in the linearized stability 
% analysis of Heimisson, Rudnicki and Lapusta (2021) for the inplane problem. Further the code
% includes relationships between slip and various fault fields that were
% omitted form the manuscript due to length and complexity. In order to
% obstain the anti-plane stability follow the instuction in Heimisson, Rudnicki and Lapusta (2021)

%If the code is used please cite: 
%Heimisson, E. R., Rudnicki, J., & Lapusta, N. (2021). Dilatancy and compaction of a rate-and-state fault in a 
%poroelastic medium: Linearized stability analysis. Journal of Geophysical Research: 
%Solid Earth, 126, e2021JB022071. https://doi.org/10.1029/2021JB022071
% 

%%
%Copyright 2021: Elias Rafn Heimisson

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
%and associated documentation files (the "Software"), to deal in the Software without restriction, 
%including without limitation the rights to use, copy, modify, merge, publish, distribute, 
%sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
%is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
%INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
%PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
%HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
%CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
clear

%% Set frictional parameters
%Independent:
a = 0.01; %rate dependence of friction
b = 0.015; % state dependence of friction
V0 = 1.0e-9; %steady state slip speed AND reference slip speed [m/s]
L = 10.0e-6; %state evolution distance [m]
f0 = 0.6; %steady state coefficient of friction
alp=0.0; %Linker-Dieterich constant
tau0 = 30.0e6; % background and initial uniform shear stress [Pa]
%Dependent: 
si0 = tau0/f0; %background and initial uniform effective normal stress [PA]

%% Shearing layer
%Independent:
epsilon = 1.0e-3; %layer half-thickness [m]
kappac = 1.0e-12; %Across layer mobility [m^2/Pa s] 
bfp = 0.44e-9; %isotropic fluid compressibility [1/Pa]
bnp = 6.0e-9; %isotropic void volume compressibility [1/Pa]
bgp = 1/(50.0e9); %isotropic solid gouge constitutent compressibility [1/Pa]
phi = 0.068; %Reference porosity/void volume
gamma = 1.7e-4; %Dilatancy coefficient (in Segall and Rice (1995) it was called epsilon)
%Dependent:
epsi = epsilon; %same as epsilon [m]
bfs = 5*bfp/9; %uniaxial fluid compressibility [1/Pa]
bns = 5*bnp/9; %uniaxial void volume compressibility [1/Pa]
bgs = 5*bgp/9; %uniaxial solid gouge constitutent compressibility [1/Pa]

beta = phi*(bfp + bnp); % Lumped compressibility for Segall and Rice (1995)
betas =phi*(bfs - bns); % Corresponding uniaxial lumped compressibility

%NOTE: bfs,bns,bgs, could be treated as dependent but that is not done in
%this study due to uncertainty in their values.

%% Poroelastic Bulk 
%Independent
% Westerly Granite values for nu, nuu, B (see Cheng (2016))
nu = 0.250; %Drained Poisson's ration
nuu = 0.331; %Undrained Poisson's ration
B = 0.810; %Skempton's coefficient
G = 30.0e9; %Shear modulus [Pa]
c = 1.0e-4; %Hydraulic diffusivity [m^2/s]
%Dependent:
alpB = 3*(nuu-nu)/(B*(1+nuu)*(1-2*nu)); %biot coefficient
kappa = c/(2*G*(1+nu)*B / (3*alpB*(1-alpB*B)*(1-2*nu) ) ); %bulk mobility m^2/(Pa s)
%% Implement equations:
% Various equations are written up using Matlab's anonymous functions

%Root finding alogrithm options:
options = optimoptions(@fsolve,'TolFun',1e-12,'TolX',1.0e-12,'MaxFunEvals',2000,'MaxIter',1000,'Display','off');


FF = @(k)kappac/(kappa*epsilon*abs(k)); %Flux non-dimensional Group, equation 18
F =  @(k) kappac/(kappa*epsilon*abs(k)); %Same
sq = @(s,k)(sqrt(1+s/(c*k^2)));  %A square root that shows up frequently
H1 = @(s,k)(1 - 2*c*k^2/s*(nuu-nu)/(1 - nu)*(1+FF(k))/(FF(k) + sq(s,k))*(sq(s,k)-1)); %Equation 16
H2 = @(s,k)((sq(s,k)-1)/(sq(s,k)+FF(k))); %Equation 17

%C3 - C7 just represent various terms that show up to shorten expressions.
%These don't have very transparent physical meanings. 
C3 = @(s,k) (F(k))/(F(k)+1)*(H2(s,k)-1);  
C4 = @(s,k) 2*epsi*abs(k)/3*B*H2(s,k)*(1+nuu)/(1-nuu);
C5 = @(s,k) kappac/(beta*epsi^2*s);
C6 = @(s,k) 2/(3*B*(1+nuu))*F(k)/(F(k)+1)*(H1(s,k)-1);
C7 = @(s,k) 2*epsi*abs(k)/(2*(1-nuu))*H1(s,k);

% Central pore-pressure:
pc = @(s,k) -(gamma*s*(2*phi + 2*G*betas*C7(s, k) - G*beta*C4(s, k) - G*bgp*C4(s, k) - 2*G*bgs*C7(s, k) + 2*G*beta*C4(s, k)*C5(s, k) + G*bgp*phi*C4(s, k) + 2*G*bgs*phi*C7(s, k) + G*bnp*phi*C4(s, k) + 2*G*bns*phi*C7(s, k) - 2))/((V0 + L*s)*(beta*phi - beta - 2*betas*C6(s, k) + beta*C3(s, k) - 2*beta*C5(s, k) - 2*beta*C3(s, k)*C5(s, k) + 2*betas*phi*C6(s, k) - beta*phi*C3(s, k) + 2*beta*phi*C5(s, k) + 2*beta*phi*C3(s, k)*C5(s, k) - G*betas*bgp*C7(s, k) - G*beta*bgs*C7(s, k) + G*betas*bgp*C3(s, k)*C7(s, k) - G*betas*bgp*C4(s, k)*C6(s, k) - 2*G*beta*bgp*C4(s, k)*C5(s, k) + G*beta*bgs*C3(s, k)*C7(s, k) - G*beta*bgs*C4(s, k)*C6(s, k) - 2*G*beta*bgs*C5(s, k)*C7(s, k) + G*betas*bgp*phi*C7(s, k) + G*beta*bgs*phi*C7(s, k) + G*betas*bnp*phi*C7(s, k) + G*beta*bns*phi*C7(s, k) - G*betas*bgp*phi*C3(s, k)*C7(s, k) + G*betas*bgp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgp*phi*C4(s, k)*C5(s, k) - G*beta*bgs*phi*C3(s, k)*C7(s, k) + G*beta*bgs*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C5(s, k)*C7(s, k) - G*betas*bnp*phi*C3(s, k)*C7(s, k) + G*betas*bnp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bnp*phi*C4(s, k)*C5(s, k) - G*beta*bns*phi*C3(s, k)*C7(s, k) + G*beta*bns*phi*C4(s, k)*C6(s, k) + 2*G*beta*bns*phi*C5(s, k)*C7(s, k) - 2*G*beta*bgs*C3(s, k)*C5(s, k)*C7(s, k) + 2*G*beta*bgs*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bgs*phi*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bns*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bns*phi*C4(s, k)*C5(s, k)*C6(s, k)));
% Average pore-pressure (used to compute the effective normal stress):
pm = @(s,k) -(gamma*s*(phi + C3(s, k) - phi*C3(s, k) + G*betas*C7(s, k) - G*bgs*C7(s, k) - G*betas*C3(s, k)*C7(s, k) + G*betas*C4(s, k)*C6(s, k) + 2*G*beta*C4(s, k)*C5(s, k) + G*bgs*C3(s, k)*C7(s, k) - G*bgs*C4(s, k)*C6(s, k) + G*bgs*phi*C7(s, k) + G*bns*phi*C7(s, k) - G*bgs*phi*C3(s, k)*C7(s, k) + G*bgs*phi*C4(s, k)*C6(s, k) - G*bns*phi*C3(s, k)*C7(s, k) + G*bns*phi*C4(s, k)*C6(s, k) - 1))/((V0 + L*s)*(beta*phi - beta - 2*betas*C6(s, k) + beta*C3(s, k) - 2*beta*C5(s, k) - 2*beta*C3(s, k)*C5(s, k) + 2*betas*phi*C6(s, k) - beta*phi*C3(s, k) + 2*beta*phi*C5(s, k) + 2*beta*phi*C3(s, k)*C5(s, k) - G*betas*bgp*C7(s, k) - G*beta*bgs*C7(s, k) + G*betas*bgp*C3(s, k)*C7(s, k) - G*betas*bgp*C4(s, k)*C6(s, k) - 2*G*beta*bgp*C4(s, k)*C5(s, k) + G*beta*bgs*C3(s, k)*C7(s, k) - G*beta*bgs*C4(s, k)*C6(s, k) - 2*G*beta*bgs*C5(s, k)*C7(s, k) + G*betas*bgp*phi*C7(s, k) + G*beta*bgs*phi*C7(s, k) + G*betas*bnp*phi*C7(s, k) + G*beta*bns*phi*C7(s, k) - G*betas*bgp*phi*C3(s, k)*C7(s, k) + G*betas*bgp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgp*phi*C4(s, k)*C5(s, k) - G*beta*bgs*phi*C3(s, k)*C7(s, k) + G*beta*bgs*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C5(s, k)*C7(s, k) - G*betas*bnp*phi*C3(s, k)*C7(s, k) + G*betas*bnp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bnp*phi*C4(s, k)*C5(s, k) - G*beta*bns*phi*C3(s, k)*C7(s, k) + G*beta*bns*phi*C4(s, k)*C6(s, k) + 2*G*beta*bns*phi*C5(s, k)*C7(s, k) - 2*G*beta*bgs*C3(s, k)*C5(s, k)*C7(s, k) + 2*G*beta*bgs*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bgs*phi*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bns*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bns*phi*C4(s, k)*C5(s, k)*C6(s, k)));
% Average of p^+ and p^- i.e. 0.5*(p^+ + p^-), introduced to simplify
% expressions:
pave = @(s,k) -(gamma*s*(2*C3(s, k) - 2*phi*C3(s, k) + G*beta*C4(s, k) + G*bgp*C4(s, k) - 2*G*betas*C3(s, k)*C7(s, k) + 2*G*betas*C4(s, k)*C6(s, k) + 2*G*beta*C4(s, k)*C5(s, k) + 2*G*bgs*C3(s, k)*C7(s, k) - 2*G*bgs*C4(s, k)*C6(s, k) - G*bgp*phi*C4(s, k) - G*bnp*phi*C4(s, k) - 2*G*bgs*phi*C3(s, k)*C7(s, k) + 2*G*bgs*phi*C4(s, k)*C6(s, k) - 2*G*bns*phi*C3(s, k)*C7(s, k) + 2*G*bns*phi*C4(s, k)*C6(s, k)))/((V0 + L*s)*(beta*phi - beta - 2*betas*C6(s, k) + beta*C3(s, k) - 2*beta*C5(s, k) - 2*beta*C3(s, k)*C5(s, k) + 2*betas*phi*C6(s, k) - beta*phi*C3(s, k) + 2*beta*phi*C5(s, k) + 2*beta*phi*C3(s, k)*C5(s, k) - G*betas*bgp*C7(s, k) - G*beta*bgs*C7(s, k) + G*betas*bgp*C3(s, k)*C7(s, k) - G*betas*bgp*C4(s, k)*C6(s, k) - 2*G*beta*bgp*C4(s, k)*C5(s, k) + G*beta*bgs*C3(s, k)*C7(s, k) - G*beta*bgs*C4(s, k)*C6(s, k) - 2*G*beta*bgs*C5(s, k)*C7(s, k) + G*betas*bgp*phi*C7(s, k) + G*beta*bgs*phi*C7(s, k) + G*betas*bnp*phi*C7(s, k) + G*beta*bns*phi*C7(s, k) - G*betas*bgp*phi*C3(s, k)*C7(s, k) + G*betas*bgp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgp*phi*C4(s, k)*C5(s, k) - G*beta*bgs*phi*C3(s, k)*C7(s, k) + G*beta*bgs*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C5(s, k)*C7(s, k) - G*betas*bnp*phi*C3(s, k)*C7(s, k) + G*betas*bnp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bnp*phi*C4(s, k)*C5(s, k) - G*beta*bns*phi*C3(s, k)*C7(s, k) + G*beta*bns*phi*C4(s, k)*C6(s, k) + 2*G*beta*bns*phi*C5(s, k)*C7(s, k) - 2*G*beta*bgs*C3(s, k)*C5(s, k)*C7(s, k) + 2*G*beta*bgs*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bgs*phi*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bns*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bns*phi*C4(s, k)*C5(s, k)*C6(s, k)));
% Changes in normal stress (not effective normal stress)
sigyy = @(s,k) (gamma*s*(2*C6(s, k) - 2*phi*C6(s, k) + G*beta*C7(s, k) + G*bgp*C7(s, k) - G*beta*C3(s, k)*C7(s, k) + G*beta*C4(s, k)*C6(s, k) + 2*G*beta*C5(s, k)*C7(s, k) - G*bgp*C3(s, k)*C7(s, k) + G*bgp*C4(s, k)*C6(s, k) - G*bgp*phi*C7(s, k) - G*bnp*phi*C7(s, k) + 2*G*beta*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*C4(s, k)*C5(s, k)*C6(s, k) + G*bgp*phi*C3(s, k)*C7(s, k) - G*bgp*phi*C4(s, k)*C6(s, k) + G*bnp*phi*C3(s, k)*C7(s, k) - G*bnp*phi*C4(s, k)*C6(s, k)))/((V0 + L*s)*(beta*phi - beta - 2*betas*C6(s, k) + beta*C3(s, k) - 2*beta*C5(s, k) - 2*beta*C3(s, k)*C5(s, k) + 2*betas*phi*C6(s, k) - beta*phi*C3(s, k) + 2*beta*phi*C5(s, k) + 2*beta*phi*C3(s, k)*C5(s, k) - G*betas*bgp*C7(s, k) - G*beta*bgs*C7(s, k) + G*betas*bgp*C3(s, k)*C7(s, k) - G*betas*bgp*C4(s, k)*C6(s, k) - 2*G*beta*bgp*C4(s, k)*C5(s, k) + G*beta*bgs*C3(s, k)*C7(s, k) - G*beta*bgs*C4(s, k)*C6(s, k) - 2*G*beta*bgs*C5(s, k)*C7(s, k) + G*betas*bgp*phi*C7(s, k) + G*beta*bgs*phi*C7(s, k) + G*betas*bnp*phi*C7(s, k) + G*beta*bns*phi*C7(s, k) - G*betas*bgp*phi*C3(s, k)*C7(s, k) + G*betas*bgp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgp*phi*C4(s, k)*C5(s, k) - G*beta*bgs*phi*C3(s, k)*C7(s, k) + G*beta*bgs*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C5(s, k)*C7(s, k) - G*betas*bnp*phi*C3(s, k)*C7(s, k) + G*betas*bnp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bnp*phi*C4(s, k)*C5(s, k) - G*beta*bns*phi*C3(s, k)*C7(s, k) + G*beta*bns*phi*C4(s, k)*C6(s, k) + 2*G*beta*bns*phi*C5(s, k)*C7(s, k) - 2*G*beta*bgs*C3(s, k)*C5(s, k)*C7(s, k) + 2*G*beta*bgs*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bgs*phi*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bns*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bns*phi*C4(s, k)*C5(s, k)*C6(s, k)));
% Fault opening and closing
dy = @(s,k) -(2*epsi*gamma*s*(beta + bgp - bgp*phi - bnp*phi + 2*betas*C6(s, k) - beta*C3(s, k) + 2*beta*C5(s, k) - bgp*C3(s, k) - 2*bgs*C6(s, k) + 2*beta*C3(s, k)*C5(s, k) + bgp*phi*C3(s, k) + 2*bgs*phi*C6(s, k) + bnp*phi*C3(s, k) + 2*bns*phi*C6(s, k)))/((V0 + L*s)*(beta*phi - beta - 2*betas*C6(s, k) + beta*C3(s, k) - 2*beta*C5(s, k) - 2*beta*C3(s, k)*C5(s, k) + 2*betas*phi*C6(s, k) - beta*phi*C3(s, k) + 2*beta*phi*C5(s, k) + 2*beta*phi*C3(s, k)*C5(s, k) - G*betas*bgp*C7(s, k) - G*beta*bgs*C7(s, k) + G*betas*bgp*C3(s, k)*C7(s, k) - G*betas*bgp*C4(s, k)*C6(s, k) - 2*G*beta*bgp*C4(s, k)*C5(s, k) + G*beta*bgs*C3(s, k)*C7(s, k) - G*beta*bgs*C4(s, k)*C6(s, k) - 2*G*beta*bgs*C5(s, k)*C7(s, k) + G*betas*bgp*phi*C7(s, k) + G*beta*bgs*phi*C7(s, k) + G*betas*bnp*phi*C7(s, k) + G*beta*bns*phi*C7(s, k) - G*betas*bgp*phi*C3(s, k)*C7(s, k) + G*betas*bgp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgp*phi*C4(s, k)*C5(s, k) - G*beta*bgs*phi*C3(s, k)*C7(s, k) + G*beta*bgs*phi*C4(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C5(s, k)*C7(s, k) - G*betas*bnp*phi*C3(s, k)*C7(s, k) + G*betas*bnp*phi*C4(s, k)*C6(s, k) + 2*G*beta*bnp*phi*C4(s, k)*C5(s, k) - G*beta*bns*phi*C3(s, k)*C7(s, k) + G*beta*bns*phi*C4(s, k)*C6(s, k) + 2*G*beta*bns*phi*C5(s, k)*C7(s, k) - 2*G*beta*bgs*C3(s, k)*C5(s, k)*C7(s, k) + 2*G*beta*bgs*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bgs*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bgs*phi*C4(s, k)*C5(s, k)*C6(s, k) + 2*G*beta*bns*phi*C3(s, k)*C5(s, k)*C7(s, k) - 2*G*beta*bns*phi*C4(s, k)*C5(s, k)*C6(s, k)));

%The characteristic equation the needs to be solved to determine stability:
CharEQ = @(s,k)( (a*si0/V0*s^2 + s*(a-b)*si0/L) + (s + V0/L)*G*abs(k)*H1(s,k)/(2*(1-nuu)) ...
      -(f0*(s + V0/L)-alp*s)*( sigyy(s,k) + pm(s,k) ) );

  
%Antiplane critical wavenumber used for initial guess in root finding algortihm.  
kcra = 2*si0*(b-a)/(L*G);
%Undrained in-plane critical wavenumber used for determining if separation of scales is honored: 
kcr = 2*si0*(b-a)*(1-nuu)/(L*G);

if kcr*epsi > 0.05
    disp('warning: separation of scales may not be honored')
    error('stop')
end

fun = @(S)([real(CharEQ(0*S(1) - 1i*S(2),S(1)));imag(CharEQ(0*S(1) - 1i*S(2),S(1)))]);

[val,fval] = fsolve(fun,[kcra V0/L],options);
kcrit = val(1) %Critical wavenumber given the parameters above
lambdacrit = 2*pi/kcrit %Corresponding critical wavelength
angular_frequency = val(2); %The angular frequency at which a critical wavelength oscillates. 

%Note, in order to obtain all solutions the user should use many
%starting points for the solver. It is possible that the solver may stop at
%a location that isn't a solution. The initial point here has been picked
%because it generally works, however be critical of the solution you obtain
%and try to validate.

