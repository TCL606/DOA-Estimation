clear all; close all; clc;

rng(1888);

c = 3e8;
M = 10;
fs = 1e9;
N = 1024;
t = linspace(0, (N - 1) / fs, N);

theta = [-10, 20, 50, 70, 80] / 180 * pi;
f = [1, 3, 5, 6, 7] * 1e6;
f0 = 500e6;
lambda = c / f0;
s = exp(2 * pi * f' * 1j * t);

snr = 12;
gm = gmdistribution(0, 1 / db2pow(snr));
gwn = zeros(M, N);
for i = 1: 1: M
    gwn(i, :) = random(gm, N)' + 1j * random(gm, N)';
end

dist = [0.25, 0.3, 0.5, 0.8, 1, 1.2] * lambda;

for dist_idx = 1: 1: length(dist)
    A = exp([0: 1: M - 1]' * 2 * pi * dist(dist_idx) / lambda * sin(theta) * -1j);
    x = A * s + gwn;
    search_points = 1024;
    theta_search = linspace(-pi / 2, pi / 2, search_points);
    a_search = exp([0: 1: M - 1]' * 2 * pi * dist(dist_idx) / lambda * sin(theta_search) * -1j);
    threhold = 0.01;

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
    subplot(length(dist), 1, dist_idx);
    plot(theta_search_degree_music, P_music);
    hold on;
    scatter(detect_theta_degree_music, P_music(peak_idx));
    title(strcat("d/\lambda=", num2str(dist(dist_idx) / lambda)));
end
xlim([-90, 90]);
