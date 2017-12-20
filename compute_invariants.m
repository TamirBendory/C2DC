% This function computes the invariants for SO(3) according to Rici's paper

% input: SH_matrix  - an array of size degree+1. The lth cell is a vector of size
% 2l+1

function [M1,M2,M3] = compute_invariants(SH_matrix)

degree = size(SH_matrix,1) - 1;

%% first order invariant (mean)
M1 = SH_matrix{1}(1);

%% second-order invariant
M2 = zeros(degree+1,1);
for i = 0:degree
    M2(i+1) =  SH_matrix{i+1}(:,1)'*SH_matrix{i+1}(:,1);
end

%% third-order invariant

% uploading the table of CGC coefficients up to degree 20
C = load('CGC_table_20'); C = C.C;

% Based on Risi's paper, compute the 4D array g
g = cell(degree+1,degree+1,degree+1);
M3 = zeros(degree+1,degree+1,degree+1);

for l1 = 0: degree
    for l2 = 0: degree
        for l = 0: degree
            g{l1+1,l2+1,l+1} = zeros(2*l+1,1);
            for m=-l:l
                for m1 = -l1:l1
                    if abs(m-m1)>l2
                        continue;
                    else
                        g{l1+1,l2+1,l+1}(m+l+1) = g{l1+1,l2+1,l+1}(m+l+1) +...
                            C{l1+1,l2+1,l+1}(l1+1+m1,l2+1+m-m1,l+1+m)*SH_matrix{l1+1}(l1+1+m1)*SH_matrix{l2+1}(l2+1+m-m1);
                    end
                end
            end 
            M3(l1+1,l2+1,l1+1) = g{l1+1,l2+1,l+1}'*SH_matrix{l+1};
        end
    end
end

