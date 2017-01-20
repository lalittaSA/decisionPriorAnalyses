function [out] = permuteBinary(n)
%permuteBinary generates all permutations of ones and zeros in N digits

z=zeros(1,n);
one_n=zeros(n,n);
for i=1:n
    for j=1:i
        one_n(i,j) = 1;
    end
end
out = z;
for i = 1:n
    temp = one_n(i,:);
    temp = perms(temp);
    one_i = unique(temp,'rows');
    out = [out; one_i];
end

out = flipud(out);
end

