% הגדרות כלליות
N_s = 1000;
theta_deg = linspace(-70, 70, N_s);  % תחום מצומצם כדי למנוע בעיות שורש
theta_rad = deg2rad(theta_deg);

% פרמטרים של התווך
x_0 = 2;
a = 1;
n_0 = 1.5;
n_x0 = sqrt(x_0/a + n_0^2);

% ביטוי לבטא
beta = n_x0 * sin(theta_rad);

% נוודא שאין שורש מרוכב ע"י פילטור
root_arg = x_0/a + n_0^2 - beta.^2;
validIdx = root_arg > 0;  % נשאר רק עם ערכים חוקיים לשורש

beta = beta(validIdx);
theta_rad = theta_rad(validIdx);  
v = -(beta.^2) ./ sqrt(root_arg(validIdx));

% חישוב נקודות הקוסטיק
x_caustic = a .* (v.^2 - n_0^2 + beta.^2);
z_caustic = 2 .* a .* beta .* (v - sqrt(n_0^2 + x_0/a - beta.^2));

% מיון לפי בטא  כדי לשמור על סדר פיזיקלי
[beta_sorted, sortIdx] = sort(beta);
x_caustic = x_caustic(sortIdx);
z_caustic = z_caustic(sortIdx);

% ציור הקוסטיק
figure; hold on;
plot(x_caustic, z_caustic, 'LineWidth', 1.4, 'Color', '#37AFE1', 'DisplayName', 'Caustic');
xlabel('x');
ylabel('z(x)');
title('Caustic Curve');
legend show;
grid on;