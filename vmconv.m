% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% vectorized upper triangle mat or convert vectorized to mat
% ALL RIGHTS RESERVED @ 2020 HAMED HONARI - JHU
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


function [M] = vmconv(v,option)

switch option 
    case 'mat2vec'
        r = size(v,1);
        % Finding the corresponding indices
        indx = nchoosek(1:r,2);
        
        % vectorizing the elements row-wise (taking elements of each rows
        % and staking them)
        for i = 1:length(indx)
            M(i) = v(indx(i,1),indx(i,2));
        end
        
    case 'vec2mat'
        % finding the size of the matrix
        r = roots([1 -1 -2*length(v)]);
        r = int32(r(r>0));

        % Finding the corresponding indices
        indx = nchoosek(1:r,2);

        
        % Now putting vectorize matrix back to a matrix

        M = zeros(r);
        for i = 1:length(indx)
            M(indx(i,1),indx(i,2)) = v(i);
        end

        M = M + M' + eye(size(M));

end
end

