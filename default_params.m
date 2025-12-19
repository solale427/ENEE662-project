function p = default_params()
% Typical EKV parameters
p.kappa = 0.7;      % gate coupling factor
p.VT0   = 0.45;   %  threshold (V) - assumed
p.UT    = 0.0259;  % thermal voltage at 300K (V)
p.Is    = 1e-6;   %  current scale (A) - assumed
p.lambda = 0.05;   % channel-length modulation - can be excludd 1/V (optional)
end
