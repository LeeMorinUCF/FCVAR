function [ db, alpha, Gamma ] = SEvec2mat( param, k, r, p )

% function param = SEmat2vec( db, alpha, Gamma, k, r, p )
% function [ db, alpha, Gamma ] = SEvec2mat( param, k, r, p )
% 
% SEmat2vec( db, alpha, Gamma, k, r, p ) and SEvec2mat( param, k, r, p ) are inverse
%       functions of each other.
% 
% SEmat2vec( db, alpha, Gamma, k, r, p ) transforms parameters in the form of
%       matrices into a vector param.
% SEvec2mat( param, k, r, p ) transforms a vector param into the
%       component matrices alpha and Gamma and the vector db.
%       
%       param is either a (1 + p*r + p^2*k) vector or a (2 + p*r + p^2*k)
%       vector depending on the size of db.
% 
% No function dependencies.
% 
%______________________________________________________



if k > 0
    Gamma = reshape( param(end - p*p*k + 1 : end), p, p*k);
    param = param(1 : end - p*p*k );
else
    Gamma = [];
end

if r > 0
    alpha = reshape( param(end - p*r + 1 : end), p, r);
    param = param(1 : end - p*r );
else
    alpha = [];
end

db = param;




% end