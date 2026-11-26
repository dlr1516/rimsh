function [m] = sh_rot_index(M, r, c)
%SH_ROT_INDEX  Entry of a band matrix addressed with centred indices.
%
%   v = SH_ROT_INDEX(M, r, c)
%
% Returns the value of the matrix M in row r and column c where indices 
% r and c have 0 value in center index.
%

    [rnum, cnum] = size(M);
    rnum = floor((rnum-1)/2); 
    cnum = floor((cnum-1)/2); 
    if (abs(r) > rnum || abs(c) > cnum)
        error(['invalid index r = ' num2str(r) '  abs(r) < ' num2str(rnum) ...
            ' or c ' num2str(c) ' abs(c) < ' num2str(cnum)])
    end
    m = M(r + rnum + 1, c + cnum + 1);
end
