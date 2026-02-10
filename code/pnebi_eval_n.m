function [pnebis] = pnebi_eval_n(n,x)

BIG_NUM = 1.0e+10;
SMALL_NUM = 1.0e-10;
%     factor, seqPrev, seqCurr, seqNext;


pnebis = zeros(n+1,1);

% if x < 0.0
%     x = -x;
% end
x = abs(x);

% If x~=0, then BesselI(0,x) = 1.0 and BesselI(k,x) = 0.0 for k > 0.
% Thus, PNEBI(0,x) = 2.0 and PNEBI(k,x) = 0.0 for k > 0.
if (x < 1e-6)
    pnebis = zeros(size(pnebis));
    pnebis(1) = 2.0;
    return;
end

% pnebis(1) = 2.0;

% Computes bessel function using back recursion
factor = 2.0 / x;
seqPrev = 0.0;  % bip
seqCurr = 1.0;  % bi
%     seqNext = 0.0;  % bim
startK = 2 * (n + floor(sqrt(40.0 * n)));
for k = startK:-1:0

    seqNext = seqPrev + factor * k * seqCurr;
    seqPrev = seqCurr;
    seqCurr = seqNext;


    if (k <= n)
        pnebis(k+1) = seqPrev;
    end

    % To avoid overflow!
    if (seqCurr > BIG_NUM)
        seqPrev = seqPrev*SMALL_NUM;
        seqCurr = seqCurr*SMALL_NUM;

        pnebis = pnebis*SMALL_NUM;

        % std::cerr << __FILE__ << "," << __LINE__ << ": ANTI-OVERFLOW!" << std::endl;

    end
end

scaleFactor = pnebi_eval_0(x) / pnebis(1);
pnebis = scaleFactor * pnebis;

end %function
