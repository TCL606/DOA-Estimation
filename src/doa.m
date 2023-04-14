clear all; close all; clc;

rng(1888);

c = 3e8;
d = 0.3;
M = 5;
fs = 1e9;
N = 1024;
t = linspace(0, (N - 1) / fs, N);

theta = [-10, 20, 50, 70, 80] / 180 * pi;
f = [1, 3, 5, 6, 7] * 1e6;
f0 = 500e6;
lambda = c / f0;
A = exp([0: 1: M - 1]' * 2 * pi * d / lambda * sin(theta) * -1j);
s = exp(2 * pi * f' * 1j * t);
