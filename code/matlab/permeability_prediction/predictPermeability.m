function p = predictPermeability(X)
%  Predict the permeability of a molecule based on mecular descriptors.
% --------------------
% Lipinski's rule of 5
% --------------------
% No more than 5 hydrogen bond donors
% No more than 10 hydrogen bond acceptors
% molecular mass less than 500 daltons
% logP that does not exceed 5
%
% Extensions:
% logP between -0.4 and 5.6
% molecular weight between 180 and 480
% 10 or fewer rotatable bonds
% TPSA less than 140 A
% --------------------
% 
% Input:
%       table X:       table containing the descriptors (columns)
% 
%                       xlogp       (octanol/water partition coefficient)
%                       mw          (molecular weight)
%                       polararea   (polar surface area in Angström)
%                       hbonddonor  (number of hydrongen bond donors)
%                       hbondacc    (number of hydroge bond acceptors)
%                       rotbonds    (number of rotatable bonds)
% Output:
%       logical p:      if p(i)=1, the metabolite has a high probability of
%                       passing a cellular membrane
          
% thresholds

t_logP_lb = -0.4;
t_logP_ub = 5.6;
t_mw = 500;
t_psa = 140;
t_hbd = 5;
t_hba = 10;
t_rot = 10;

p = ...
    X.xlogp >= t_logP_lb & ...
    X.xlogp <= t_logP_ub & ...
    X.mw <= t_mw & ...
    X.polararea <= t_psa & ...
    X.hbonddonor <= t_hbd & ...
    X.hbondacc <= t_hba & ...
    X.rotbonds <= t_rot;

end





