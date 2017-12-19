% WIGNER_D Calculate a Wigner D-matrix
%
% Usage
%    D = wigner_d(j, angles);
%
% Input
%    j: A positive integer corresponding to the order of the matrix.
%    angles: A 3-by-n array consisting of the Euler angles for which we want
%       to calculate the Wigner D-matrices. These are specified in the z-y-z
%       convention as active rotations.
%
% Output
%    D: A (2j+1)-by-(2j+1)-by-n array consisting of the Wigner D-matrices of
%       order j with the specified Euler angles. The rows correspond to the
%       m' indices while the columns correspond to m indices. Both run from
%       -j to j.
%
% Note
%    The easyspin package uses different conventions for its 'wignerd'
%    function. Notably, the order of the indices run from j to -j and the
%    matrix is multiplied by (-1)^j. To get the same matrices, you would
%    therefore calculate
%
%       (-1)^j*flip(flip(wigner_d(j, angles), 1), 2) .

% Author: Joakim Anden, janden@math.princeton.edu

function D = wigner_d(j, angles)
	% Calculate Wigner's (small) d-matrix first from beta.
	D = small_d(j, angles(2,:));

	% Form Wigner's D-matrix by multiplying in the alpha and gamma angles.
	D = bsxfun(@times, D, permute(exp(-i*[-j:j]'*angles(1,:)), [1 3 2]));
	D = bsxfun(@times, D, permute(exp(-i*[-j:j]'*angles(3,:)), [3 1 2]));
end

function d = small_d(j, beta)
	% Make sure beta is a column vector.
	beta = beta(:);

	n = size(beta, 1);

	% We only ever need powers of these cosines and sines, so precompute them
	% and their squares.
	co = cos(beta/2);
	si = sin(beta/2);

	cosq = co.^2;
	sisq = si.^2;

	% Initialize the output d-matrix. The first two indices are mp and m,
	% while the third is the index of the beta angle.
	d = zeros(2*j+1, 2*j+1, n);

	% To speed up calculations, we keep track of two sets of coefficients
	% a_{mp,m}, b_{mp,m} as we loop through m and mp. The first are given
	% by
	%
	%    a_{mp,m} = [(j+mp+1)*...*(j-mp)]*[(j-m+1)*...*(j+m)] ,
	%
	% while the second are given by
	%
	%    b_{mp,m} = 1/(m-mp)! .
	%
	% During the innermost (t) loop, a third set of coefficients is used
	% whose formula is given by
	%
	%   c_{mp,m,t} = a_{mp,m}^{0.5} * b_{mp,m} (-1)^t * [(j+mp-t+1)*...*(j+mp)]
	%      * [(j-m-t+1)*...*(j-m)] / [(m-mp+1)*...*(m-mp+t)] / t! .

	% Initialize our coefficients that will be updated during the loops. The
	% zeroth level (a0, b0, co0, si0) are updated during the m loop and used to
	% initialize the first level (a1, b1, co1, si1), which is updated during
	% the mp loop, and so on.
	a0 = 1;
	b0 = 1;
	co0 = co.^(2*j);
	si0 = 1;

	for m = -j:j
		a1 = a0;
		b1 = b0;
		co1 = co0;
		si1 = si0;

		% Due to symmetries, we only need to calculate one triangle of the
		% d-matrix explicitly, so we restrict mp < -abs(m).
		for mp = -j:-abs(m)
			% For a fixed pair (m, mp), the entry in the d-matrix is given by a
			% linear combination of sine and cosine powers. The coefficients are
			% given by c_{m,mp,t}. We initialize this for t = 0 using the values
			% of a_{m,mp} and b_{m,mp} and update at each iteration of the t
			% loop.
			c2 = sqrt(a1)*b1;
			co2 = co1;
			si2 = si1;

			val = zeros(n, 1);
			for t = 0:j+mp
				val = val + c2*(co2.*si2);

				c2 = -c2*(j+mp-t)*(j-m-t)/((m-mp+t+1)*(t+1));
				co2 = co2./cosq;
				si2 = si2.*sisq;
			end

			d(j+mp+1,j+m+1,:) = permute(val, [2 3 1]);

			% Update first-level coefficients.
			a1 = a1/((j+mp+1)*(j-mp));
			b1 = b1*(m-mp);
			co1 = co1.*co;
			si1 = si1./si;
		end

		% Update zeroth-level coefficients.
		a0 = a0*(j-m)*(j+m+1);
		b0 = b0/(m+1+j);
		co0 = co0./co;
		si0 = si0.*si;
	end

	% The symmetries of the d-matrix lets us fill in the rest by first flipping
	% across the antidiagonal.
	temp = flip(d, 2);
	temp = reshape(temp, [(2*j+1)^2 n]);
	temp(1:(2*j+1)+1:end,:) = 0;
	temp = reshape(temp, [(2*j+1)*ones(1, 2) n]);
	temp = flip(permute(temp, [2 1 3]), 2);
	d = d+temp;

	% The second symmetry is across the diagonal, with sign changes according to
	% the parity of m and mp.
	temp = permute(d, [2 1 3]);
	temp = reshape(temp, [(2*j+1)^2 n]);
	temp(1:(2*j+1)+1:end,:) = 0;
	temp = reshape(temp, [(2*j+1)*ones(1, 2) n]);
	temp = bsxfun(@times, temp, (-1).^(0:2*j)'*(-1).^(0:2*j));
	d = d+temp;
end

