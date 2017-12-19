% CLEBSCH_GORDAN Calculate a Clebsch-Gordan coefficient
%
% Usage
%    c = clebsch_gordan(j1, j2, j, m1, m2, m);
%
% Input
%    j1, j2, j, m1, m2, m: The indices of the coefficient. These must satisfy
%       j1, j2, j greater than or equal to zero, and m1, m2, m smaller than or
%       equal to j1, j2, and j, respectively.
%
% Output
%    c: The Clebsch-Gordan coefficient <j1,m1,j2,m2|j,m>.

% Author: Joakim Anden, janden@math.princeton.edu

function c = clebsch_gordan(j1, j2, j, m1, m2, m)
	% If m is not equal to m1+m2, the coefficient is always zero. Likewise, j
	% needs to be between abs(j1-j2) and j1+j2.
	if m ~= m1+m2 || j < abs(j1-j2) || j > j1+j2
		c = 0;
		return;
	end

	% To simplify calculations, we want to assume that j1 >= j2, so otherwise
	% we flip and use the following relation.
	if j1 < j2
		c = (-1)^(j-j1-j2)*clebsch_gordan(j2, j1, j, m2, m1, m);
		return;
	end

	% To calculate the Clebsch-Gordan coefficients for (j1, j2, j, m1, m2, m)
	% we first calculate them for a number of parameters (j1, j2, j, n1, n2, n)
	% and use recurrence formulas to arrive at the correct coefficient.

	% First, we calculate the coefficients for n = j. These are given for n1
	% ranging from j-j2 to j1 and n2 ranging from j2 to j-j1. We set the last
	% coefficient, corresponding to (n1 = j1, n2 = j-j1) to 1, use a
	% recurrence formula for the rest, then normalize all coefficients to have
	% squares sum to one.
	cj = zeros(1, j2+j1-j+1);
	cj(end) = 1;

	% Start at the end and work our way backwards.
	k = numel(cj);
	n2 = j-j1;
	n1 = j1;

	while n2 < j2
		% Use recurrence formula to calculate next coefficient.
		cj(k-1) = -sqrt(((j2-n2)*(j2+n2+1))/((j1-n1+1)*(j1+n1)))*cj(k);
		n1 = n1-1;
		n2 = n2+1;
		k = k-1;
	end

	% Normalize to have sum of squares one.
	cj = cj/sqrt(sum(abs(cj).^2));

	% In cj, n1 runs from j-j2 to j1, but we need those from m1 to m1+(j-m).
	% We therefore cut off cj or add zeros where necessary.
	if m1 > j-j2
		cj = cj(1+(m1-(j-j2)):end);
	else
		cj = [zeros(1, (j-j2)-m1) cj];
	end

	if m1+(j-m) < j1
		cj = cj(1:end-(j1-(m1+(j-m))));
	else
		cj = [cj zeros(1, (m1+(j-m))-j1)];
	end

	% Starting at n = j-1, we move downwards until n = m.
	n = j-1;
	while n >= m
		% Depending on n, certain Clebsch-Gordan coefficients will be zero
		% (if n1 < n-j2) or not useful for our target coefficient (if
		% n1 < m1), so we set a minimum for n1 accordingly.
		n1_min = n-j2;
		if n1_min < m1
			n1_min = m1;
		end

		% The same reasoning gives us a maximum for n1.
		n1_max = m1+(n-m);
		if n1_max > j1
			n1_max = j1;
		end

		% Here k is the index in the cj vector that we want to update. Since
		% we only care about a certain set of n1, we only update these values.
		for k = n1_min-m1+1:n1_max-m1+1
			n1 = m1+(k-1);
			n2 = n-n1;

			% The recurrence formula gets us the coefficient (n1, n2, n) from
			% (n1+1, n2, n+1) and (n1, n2+1, n).
			cj(k) = sqrt((j1+(n1+1))*(j1-n1))*cj(k+1) ...
				+sqrt((j2+(n2+1))*(j2-n2))*cj(k);
		end

		% Since this factor is the same for all recurrences, we calculate it
		% separately.
		cj = cj/sqrt((j+n+1)*(j-(n+1)+1));

		n = n-1;
	end

	% At the end, we have the coefficient (m1, m2, m) at the first entry.
	c = cj(1);
end

