% This script generates a table of CGC coefficients up to l_max and saves
% it

l_max = 20;
C = cell(l_max+1,l_max+1,l_max+1);

for l1 = 0:l_max
    for l2 = 0:l_max
        for l = 0:l_max
            C{l1+1,l2+1,l+1} = zeros(2*l1+1,2*l2+1,2*l+1);
            for m1 = -l1:l1
                for m2 = -l2:l2
                    for m = -l:l
                        C{l1+1,l2+1,l+1}(l1+1+m1,l2+1+m2,l+1+m) = clebsch_gordan(l1, l2, l, m1, m2, m);
                    end
                end
            end
        end
    end
end

str = strcat('CGC_table_',num2str(l_max));
save(str,'C');