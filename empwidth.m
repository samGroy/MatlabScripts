%calculates channel width for each node
q=cread('dip50_13.q',51);
secperyear=31556926;
w=10*(q/secperyear).^.5;
ctrisurf('dip50_13',51,w);
colorbar;
fprintf('minimum channel width: %.5f\n', (min(w)))
fprintf('maximum channel width: %.5f\n', (max(w)))
