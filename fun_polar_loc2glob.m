%% Inputs
% RGL : origin of the local PCS in the global PCS
% AGL : origin of the local PCS in the global PCS
% RLI : point of interest in local PCS
% ALI : point of interest in local PCS

%% Outputs
% RGI : point of interest in global PCS
% AGI : point of interest in global PCS

%% Example
% RGL = 2.8;
% AGL = pi-pi/6;
% RGI = 1.2;
% AGI = pi/6;
% RLI = 3.555277766926234;
% ALI = 2.617993877991494;

function [RGI,AGI]=fun_polar_loc2glob(RGL,AGL,RLI,ALI)

% syms RGL AGL RGI AGI ALI RLI
% e1 = RGI.*cos(AGI) == RLI.*cos(ALI)+RGL.*cos(AGL);
% e2 = RGI.*sin(AGI) == RLI.*sin(ALI)+RGL.*sin(AGL);
% e  = [e1,e2];
% s  = solve(e,RGI,AGI);
% simplify(s.RGI)
% simplify(s.AGI)

RGI = [  (cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2), ...
 -(cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2)];

AGI = [  2.*atan((RGL + RLI - RGL.*(cos(AGL) + 1) - RLI.*(cos(ALI) + 1) + (cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2))./(RGL.*sin(AGL) + RLI.*sin(ALI))), ...
 -2.*atan((RGL.*(cos(AGL) + 1) - RLI - RGL + RLI.*(cos(ALI) + 1) + (cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2))./(RGL.*sin(AGL) + RLI.*sin(ALI)))];
 
[RGI,ind_RGI] = max(RGI');
RGI           = RGI';
AGI2          = AGI(:,1);

for i0 = 1:length(AGI2)
    AGI2(i0,1)  = AGI(i0,ind_RGI(i0));
    
    if AGI2(i0,1) <0
        AGI2(i0,1) = AGI2(i0,1) + 2*pi;
    end
    
    if AGI2(i0,1)>2*pi
        AGI2(i0,1) = 2*pi-AGI2(i0,1);
    end
end

AGI = AGI2;

end
