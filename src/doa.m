clear all; close all; clc;

rng(1888);

c = 3e8;
dist = 0.3;
M = 10;
fs = 1e9;
N = 1024;
t = linspace(0, (N - 1) / fs, N);

theta = [-10, 20, 50, 70, 80] / 180 * pi;
f = [1, 3, 5, 6, 7] * 1e6;
f0 = 500e6;
lambda = c / f0;
A = exp([0: 1: M - 1]' * 2 * pi * dist / lambda * sin(theta) * -1j);
s = exp(2 * pi * f' * 1j * t);

snr = 12;
gm = gmdistribution(0, 1 / db2pow(snr));
gwn = zeros(M, N);
for i = 1: 1: M
    gwn(i, :) = random(gm, N)' + 1j * random(gm, N)';
end

x = A * s + gwn;
search_points = 1024;
theta_search = linspace(-pi / 2, pi / 2, search_points);
a_search = exp([0: 1: M - 1]' * 2 * pi * dist / lambda * sin(theta_search) * -1j);
threhold = 0.01;

% Capon
Rxx = (x * x') / N;
inv_Rxx = inv(Rxx);
P_capon = zeros(1, search_points);
for i = 1: 1: search_points
    P_capon(i) = abs(1 ./ (a_search(:, i)' * inv_Rxx * a_search(:, i)));
end
P_capon = P_capon / max(P_capon);
peak_idx = FindPeak(P_capon, 10, 0);
ignore_idx = [];
n = 0;
for i = peak_idx
    n = n + 1;
    if P_capon(i) < threhold
        ignore_idx = [n, ignore_idx];
    end
end
peak_idx(ignore_idx)= [];
theta_search_degree_capon = theta_search / pi * 180;
detect_theta_degree_capon = theta_search_degree_capon(peak_idx);
subplot(3, 1, 1);
plot(theta_search_degree_capon, P_capon);
hold on;
scatter(detect_theta_degree_capon, P_capon(peak_idx));
xlim([-90, 90]);

% MUSIC
Rxx = (x * x') / N;
[V,D] = eig(Rxx);                 
d = diag(D);
[d_sort, idx_sort] = sort(d);
noise_eig_num = 2;
idx_noise = idx_sort(1: 1 + noise_eig_num);
Vn = V(:, idx_noise);
P_music = zeros(1, search_points);
for i = 1: 1: search_points
    P_music(i) = abs(1 ./ (a_search(:, i)' * Vn * Vn' * a_search(:, i)));
end
P_music = P_music / max(P_music);
peak_idx = FindPeak(P_music, 10, 0);
ignore_idx = [];
n = 0;
for i = peak_idx
    n = n + 1;
    if P_music(i) < threhold
        ignore_idx = [n, ignore_idx];
    end
end
peak_idx(ignore_idx)= [];
theta_search_degree_music = theta_search / pi * 180;
detect_theta_degree_music = theta_search_degree_music(peak_idx);
subplot(3, 1, 2);
plot(theta_search_degree_music, P_music);
hold on;
scatter(detect_theta_degree_music, P_music(peak_idx));
xlim([-90, 90]);

% ESPRIT
[U, D, V] = svd(Rxx);
Us = U(:, 1: 5);
Usx = Us(1: M - 1,:); 
Usy = Us(2: M,:);
Psi = Usx \ Usy;
[~, phi] = eig(Psi);
phi = diag(phi);
omega = angle(phi);
detect_theta_degree_esprit = -asin(omega / 2 / pi / dist * lambda) / pi * 180;
subplot(3, 1, 3);
plot_esprit_x = cos(detect_theta_degree_esprit / 180 * pi);
plot_esprit_y = sin(detect_theta_degree_esprit / 180 * pi);
scatter(plot_esprit_x, plot_esprit_y);
xlim([-1, 1]);