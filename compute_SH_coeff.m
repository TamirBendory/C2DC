%This function computes the SH expansion of a function F, already defined on
%the sphere
% output: An array of spherical hatmonics. For each l, the array contains
% 2l+1 coefficients

function SH_matrix = compute_SH_coeff(F,degree)

k = 1; SH_matrix = cell(degree+1,1);
for l = 0:degree
    coeffs = zeros(2*l+1,1);
    for m = -l:l
        coeffs(m+l+1) = sum2(F.*spherefun.sphharm(l,m));        
        k = k + 1;
    end
    SH_matrix{l+1} = coeffs;
end
