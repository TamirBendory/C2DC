clear; close all; clc;

L = 30;
theta = linspace(0,pi);
phi = linspace(0,2*pi-1/L);
[Theta,Phi] = meshgrid(theta,phi);

y = sph_harmonic(1, 0, Theta, Phi);