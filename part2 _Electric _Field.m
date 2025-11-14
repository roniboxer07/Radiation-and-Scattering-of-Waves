clc; 
clear;

% --- קבועים ---
mu0 = 4 * pi * 1e-7;          % חדירות מגנטית של ריק [H/m]
eps0 = 8.854e-12;             % קבוע דיאלקטרי של ריק [F/m]
freq = 15e9;                  % תדר הגל [Hz]
w = 2 * pi * freq;            % תדר זוויתי [rad/s]
k0 = w * sqrt(mu0 * eps0);    % מספר גל חופשי במרחב [rad/m]
lambda_0 = 2 * pi / k0;       % אורך גל חופשי [m]
eta = sqrt(mu0 / eps0);       % אימפדנס של ריק [Ohm]
E0 = 1;                       % עוצמת שדה נומינלית
n0 = 1.5;                     % מקדם שבירה בסיסי
a = 1;                        % פרמטר של פרופיל מקדם השבירה
I0 = 1;                       % זרם המקור (מקור נקודתי)
x0 = 12;                      % מיקום יצאת הקרן במרחב (בציר x)


% --- טווחים ---
N = 1000;
theta_field_deg = linspace(-90, 90, N);
theta_field = deg2rad(theta_field_deg);
beta = n0 .* sin(theta_field).';  % עמודה

v_lin = linspace(-20, 20, N).';  % עמודה
%v = sqrt(x./a + n0^2 - beta.^2);

% --- רשת דו מימדית (v,beta) ---
[VV, BB] = meshgrid(v_lin, beta); %יוצרת מטריצות בגודל N X N  


% חישוב v0 ו-vt
v0 = sqrt(x0/a + n0^2 - BB.^2);
vt = sqrt((a * (BB.^2 - n0^2))/a + n0^2 - BB.^2);  % v_t לפי x_t
% פונקציית פאזה
phase = a * ((2/3).*vt.^3 + 2.*BB.^2 .*vt - (2/3).*v0.^3 - 2.*BB.^2 .* v0) ...
      - a * ((2/3).*VV.^3 + 2.*BB.^2 .* VV - (2/3).*vt.^3 - 2.*BB.^2 .* vt) .* (VV > vt);

% מיפוי למרחב (x,z)
x = a .* (VV.^2 - n0^2 + BB.^2);
z = 2 .* a .* BB .* (VV - v0);

% מקדם שבירה n(v,beta)
n_vb = sqrt(VV.^2 - BB.^2);

%  Ey
h = 8 .* pi .*k0.* n_vb .*4.*a.^2 .*(VV.^2+ VV.*((BB.^2 ./ v0)- v0)-BB.^2);

Ey = 1j .* w .* mu0 .* I0 .* exp(1j .* pi ./ 4) .* ...
     exp(1j .* k0 .* phase) .* sqrt(1./h);

% --- יצירת גריד להצגה ---
x_axis = linspace(-50, 250, 1000);
z_axis = linspace(-150, 150, 1000);
[x_grid, z_grid] = meshgrid(x_axis, z_axis);

% אינטרפולציה
Ey_grid = griddata(x(:), z(:), Ey(:), x_grid, z_grid); %הופך את השדה (שמחושב לפי קרניים) למפה דו-ממדית על הרשת הקבועה
Ey_grid(isnan(Ey_grid)) = 0; %NaN מתוקן ל־0 למניעת רעש

% --- גרף ---
f = figure;
surface(x_grid, z_grid, log10(abs(Ey_grid)));%log10(abs(...)) – פותח את טווח הדינאמי ומאפשר לראות מבנים חלשים בשדה.
shading interp;
colormap jet;
c = colorbar;
c.Label.String = 'log10(abs(Ey_grid))) scale';
xlabel('X');
ylabel('z(x)');
title('Electric Field');
view(2);