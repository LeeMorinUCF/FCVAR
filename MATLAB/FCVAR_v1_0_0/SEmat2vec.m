function param = SEmat2vec( db, alpha, Gamma, k, r, p )

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



param = db;

if r > 0
    param = [  param reshape( alpha, 1, p*r ) ];
end

if k > 0
    param = [  param reshape( Gamma, 1, p*p*k ) ];
end



% end