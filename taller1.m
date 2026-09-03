clc; clear; close all;

% --- Parámetros y secuencias iniciales ---
N = 8;
n = 0:N-1;
k = 0:N-1;

x = [1, 1, 0, 0, 1, 0, 0, 0];
y = [1, 0, 1, 1, 0, 0, 1, 0];

% =========================================================================
% 1. IMPLEMENTACIÓN DE LA DTFS (Análisis y Síntesis)
% =========================================================================

% Fórmulas de Análisis (Cálculo de los coeficientes X[k] y Y[k])
X = zeros(1, N); Y = zeros(1, N);
for kv = 0:N-1
    X(kv+1) = (1/N) * sum(x .* exp(-1j * (2*pi/N) * kv .* n));
    Y(kv+1) = (1/N) * sum(y .* exp(-1j * (2*pi/N) * kv .* n));
end

% Fórmulas de Síntesis (Reconstrucción de x[n] y y[n])
x_rec = zeros(1, N); y_rec = zeros(1, N);
for nv = 0:N-1
    x_rec(nv+1) = sum(X .* exp(1j * (2*pi/N) * k .* nv));
    y_rec(nv+1) = sum(Y .* exp(1j * (2*pi/N) * k .* nv));
end

% Se toma la parte real para limpiar errores numéricos de punto flotante
x_rec = real(x_rec); 
y_rec = real(y_rec);

% =========================================================================
% 2. VERIFICACIÓN NUMÉRICA DE PROPIEDADES
% =========================================================================
disp('--- Verificación de Propiedades de la DTFS ---');
disp('(Nota: Los errores numéricos del orden de 1e-15 se consideran cero)');

% A. Linealidad: DTFS(ax + by) = a*DTFS(x) + b*DTFS(y)
a = 2; b = 3; % Constantes arbitrarias
sec_lin = a*x + b*y;
DTFS_lin = zeros(1, N);
for kv = 0:N-1
    DTFS_lin(kv+1) = (1/N) * sum(sec_lin .* exp(-1j * (2*pi/N) * kv .* n));
end
DTFS_lin_prop = a*X + b*Y;
err_linealidad = max(abs(DTFS_lin - DTFS_lin_prop));
fprintf('Error de Linealidad: %e\n', err_linealidad);

% B. Desplazamiento
n0 = 2; % Desplazamiento arbitrario
x_shift = circshift(x, n0);
X_shift = zeros(1, N);
for kv = 0:N-1
    X_shift(kv+1) = (1/N) * sum(x_shift .* exp(-1j * (2*pi/N) * kv .* n));
end
% Verificación de magnitud
err_mag = max(abs(abs(X_shift) - abs(X)));
% Verificación de fase (comparando la forma exponencial compleja)
X_shift_prop = X .* exp(-1j * (2*pi/N) * k * n0);
err_fase = max(abs(X_shift - X_shift_prop));

fprintf('Error Desplazamiento (Magnitud): %e\n', err_mag);
fprintf('Error Desplazamiento (Fase/Complejo): %e\n', err_fase);

% C. Multiplicación (Convolución circular en frecuencia)
sec_mult = x .* y;
DTFS_mult = zeros(1, N);
for kv = 0:N-1
    DTFS_mult(kv+1) = (1/N) * sum(sec_mult .* exp(-1j * (2*pi/N) * kv .* n));
end
% Convolución circular usando la función nativa cconv de MATLAB
DTFS_mult_prop = cconv(X, Y, N);
err_mult = max(abs(DTFS_mult - DTFS_mult_prop));
fprintf('Error Multiplicación (Conv. Circular): %e\n', err_mult);

% Verificación de Reconstrucción Exacta
err_rec_x = max(abs(x - x_rec));
err_rec_y = max(abs(y - y_rec));
fprintf('Error Reconstrucción x[n]: %e\n', err_rec_x);
fprintf('Error Reconstrucción y[n]: %e\n\n', err_rec_y);

% =========================================================================
% 3. GRÁFICOS DE MAGNITUD, FASE Y RECONSTRUCCIÓN
% =========================================================================

% Limpieza de ruido numérico en la fase para que las gráficas sean legibles
X_phase = angle(X); X_phase(abs(X) < 1e-10) = 0;
Y_phase = angle(Y); Y_phase(abs(Y) < 1e-10) = 0;

% Figura 1: Análisis de x[n]
figure('Name', 'Análisis de x[n]', 'Color', 'w');
subplot(2,2,1); stem(n, x, 'filled', 'b'); 
title('Original x[n]'); xlabel('n'); ylim([0 1.5]); grid on;
subplot(2,2,2); stem(n, x_rec, 'filled', 'r'); 
title('Reconstruida x[n]'); xlabel('n'); ylim([0 1.5]); grid on;
subplot(2,2,3); stem(k, abs(X), 'filled', 'k'); 
title('Espectro Magnitud |X[k]|'); xlabel('k'); grid on;
subplot(2,2,4); stem(k, X_phase, 'filled', 'm'); 
title('Espectro Fase \angle X[k]'); xlabel('k'); grid on;

% Figura 2: Análisis de y[n]
figure('Name', 'Análisis de y[n]', 'Color', 'w');
subplot(2,2,1); stem(n, y, 'filled', 'b'); 
title('Original y[n]'); xlabel('n'); ylim([0 1.5]); grid on;
subplot(2,2,2); stem(n, y_rec, 'filled', 'r'); 
title('Reconstruida y[n]'); xlabel('n'); ylim([0 1.5]); grid on;
subplot(2,2,3); stem(k, abs(Y), 'filled', 'k'); 
title('Espectro Magnitud |Y[k]|'); xlabel('k'); grid on;
subplot(2,2,4); stem(k, Y_phase, 'filled', 'm'); 
title('Espectro Fase \angle Y[k]'); xlabel('k'); grid on;