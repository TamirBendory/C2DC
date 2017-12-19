% SPH_HARMONIC Calculate a spherical harmonic function
%
% Usage
%    y = sph_harmonic(j, m, theta, phi);
%
% Input
%    j, m: The indices of the spherical harmonic such that j is non-
%       negative and |m| is less than or equal to j.
%    theta, phi: The latitude and longitude angles, respectively, at which
%       the spherical harmonic is to be evaluated. These are in the form of
%       an array.
%
% Output
%    y: The value of the spherical harmonic function of index (j, m) at the
%       latitudes and longitudes specified.

% Author: Joakim Anden, janden@math.princeton.edu

function y = sph_harmonic(j, m, theta, phi)
	y = assoc_legendre(j, m, cos(theta)).*exp(i*m*phi);
	if m < 0
		y = sqrt((2*j+1)/(4*pi)*prod((j+m+1):(j-m)))*y;
	else
		y = sqrt((2*j+1)/(4*pi)/prod((j-m+1):(j+m)))*y;
	end
end

