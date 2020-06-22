%% sigma cartesian local
sigma.cartesian.local.xx  = @(x,y) (-y.*(+3*x.^2 + y.^2))./((x.^2 + y.^2).^2);
sigma.cartesian.local.yy  = @(x,y) (+y.*(+1*x.^2 - y.^2))./((x.^2 + y.^2).^2);
sigma.cartesian.local.xy  = @(x,y) (+x.*(+1*x.^2 - y.^2))./((x.^2 + y.^2).^2);

%% sigma polar local
sigma.polar.local.xx      = @(r,a) -(3*sin(a) - 2*sin(a).^3)./r;
sigma.polar.local.yy      = @(r,a) +(1*sin(a) - 2*sin(a).^3)./r;
sigma.polar.local.xy      = @(r,a) -(1*cos(a) - 2*cos(a).^3)./r;

%% sigma cartesian global
sigma.cartesian.global.xx = @(x,y,t) ...
    1*sigma.cartesian.local.xx(x,y).*cos(t).^2 + ...
    1*sigma.cartesian.local.yy(x,y).*sin(t).^2 + ...
    2*sigma.cartesian.local.xy(x,y).*sin(t).*cos(t);

sigma.cartesian.global.yy = @(x,y,t) ...
    1*sigma.cartesian.local.xx(x,y).*sin(t).^2 + ...
    1*sigma.cartesian.local.yy(x,y).*cos(t).^2 - ...
    2*sigma.cartesian.local.xy(x,y).*sin(t).*cos(t);

sigma.cartesian.global.xy = @(x,y,t) ...
    (sigma.cartesian.local.yy(x,y)-sigma.cartesian.local.xx(x,y)).*sin(t).*cos(t) + ...
    (sigma.cartesian.local.xy(x,y)).*(cos(t).^2 - sin(t).^2);

%% sigma polar global 
sigma.polar.global.xx = @(r,a,t) ...
    1*sigma.polar.local.xx(r,a).*cos(t).^2 + ...
    1*sigma.polar.local.yy(r,a).*sin(t).^2 + ...
    2*sigma.polar.local.xy(r,a).*sin(t).*cos(t);

sigma.polar.global.yy = @(r,a,t) ...
    1*sigma.polar.local.xx(r,a).*sin(t).^2 + ...
    1*sigma.polar.local.yy(r,a).*cos(t).^2 - ...
    2*sigma.polar.local.xy(r,a).*sin(t).*cos(t);

sigma.polar.global.xy = @(r,a,t) ...
    (sigma.polar.local.yy(r,a)-sigma.polar.local.xx(r,a)).*sin(t).*cos(t) + ...
    (sigma.polar.local.xy(r,a)).*(cos(t).^2 - sin(t).^2);
