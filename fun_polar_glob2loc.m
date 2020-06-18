%% Inputs
% RGL : origin of the local PCS in the global PCS
% AGL : origin of the local PCS in the global PCS
% RGI : point of interest in global PCS
% AGI : point of interest in global PCS

%% Outputs
% RLI : point of interest in local PCS
% ALI : point of interest in local PCS

%% Example
% RGL = 2.8;
% AGL = pi-pi/6;
% RGI = 1.2;
% AGI = pi/6;
% RLI = 3.555277766926234;
% ALI = 2.617993877991494;

function [RLI,ALI]=fun_polar_glob2loc(RGL,AGL,RGI,AGI,Cont_val)

if nargin<5
    Cont_val = true;
end

% syms RGL AGL RGI AGI ALI RLI
% e1 = RGI.*cos(AGI) == RLI.*cos(ALI)+RGL.*cos(AGL);
% e2 = RGI.*sin(AGI) == RLI.*sin(ALI)+RGL.*sin(AGL);
% e  = [e1,e2];
% s  = solve(e,RLI,ALI);
% simplify(s.RLI)
% simplify(s.ALI)

RLI = [-(cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2), ... 
  (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2)];

ALI = [ -2.*atan((RGL - RGI + RGI.*(cos(AGI) + 1) - RGL.*(cos(AGL) + 1) + (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2))./(RGI.*sin(AGI) - RGL.*sin(AGL))), ...
  2.*atan((RGI - RGL - RGI.*(cos(AGI) + 1) + RGL.*(cos(AGL) + 1) + (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2))./(RGI.*sin(AGI) - RGL.*sin(AGL)))];
 
% Phase control 
if Cont_val
    [RLI,ind_RLI] = max(RLI');
    RLI           = RLI';
    ALI2          = ALI(:,1);
    
    for i0 = 1:length(ALI2)
        ALI2(i0,1)  = ALI(i0,ind_RLI(i0));
        
        if ALI2(i0,1) <0
            ALI2(i0,1) = ALI2(i0,1) + 2*pi;
        end
        
        if ALI2(i0,1)>2*pi
            ALI2(i0,1) = 2*pi-ALI2(i0,1);
        end
    end
    
    ALI = ALI2;
    
end

end
