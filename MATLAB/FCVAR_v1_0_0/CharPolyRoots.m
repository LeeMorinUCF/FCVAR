function charPolyRoots = CharPolyRoots(b, alpha, beta, Gamma, print2screen)

% function charPolyRoots = CharPolyRoots(b, alpha, beta, Gamma, print2screen)
% Lee Morin and Morten Nielsen
% August 21, 2011
%
% CharPolyRoots calculates the roots of the characteristic polynomial
%       for the usual cointegrated VAR model and plots them
%       with the unit circle transformed for the fractional model
% 
% input = scalar b is the fractional differencing parameter.
%       matrices alpha, beta and Gamma are the parameters for the
%           fractional cointegrated VAR model.
%       matrix Gamma is a p x p*k matrix of coefficients in the lag
%           polynomial in delta d * Lbk * x.
%       print2screen is an indicator to determine if a graph of the roots
%           will be produced (to complement other statistics printed to screen).
% 
% output = complex vector charPolyRoots with the roots of the characteristic polynomial.
% 
% No dependencies.
% 
% Note: This function calculates roots by constructing a matrix for a 1st
%       order lag polynomial on a big vector of stacked lagged vectors.
%       The characteristic roots are then the eigenvalues of this matrix.
% 
%______________________________________________________


% Get dimension of system.
p = size(alpha, 1);
% Get lag length from Gamma.
k = floor(size(Gamma,2)/p + 0.1);


% First construct a matrix for the stacked version of the VAR 
%   with a 1st order lag polynomial.
PiStar = eye(p) + alpha*beta';

if k > 0
    Gamma1 = Gamma(:, 1 : p );
    PiStar = PiStar + Gamma1;
end
for i = 2:k
    
    Gammai = Gamma(:, (i-1)*p + 1 : i*p );
    GammaiMinus1 = Gamma(:, (i-2)*p + 1 : (i-1)*p );
    
    PiStar = [ PiStar (Gammai - GammaiMinus1) ];
    
end
if k > 0
    Gammak = Gamma(:, (k-1)*p + 1 : k*p );
    PiStar = [ PiStar  ( - Gammak ) ];
end

% Pad with an identity for the transition of the lagged variables.
PiStar = [  PiStar;...
            eye(p*k) zeros( p*k, p*(k>0) )  ];

        
% The roots are then the eigenvalues of the matrix PiStar.
charPolyRoots = eig(PiStar);


% Generate graphs if estimation results are printed to screen.
if print2screen
    
    % Now calculate the line for the transformed unit circle.
    % First do the negative half.
    unitCircle = ( pi : - 0.001 : 0)';
    psi = - (pi - unitCircle)/2;
    unitCircleX = cos( - unitCircle);
    unitCircleY = sin( - unitCircle);
    transformedUnitCircleX = [ (1 - (2*cos(psi)).^b.*cos(b*psi)) ];
    transformedUnitCircleY = [ (    (2*cos(psi)).^b.*sin(b*psi)) ];
    % Then do the positive half.
    unitCircle = (0 : 0.001 : pi)';
    psi = (pi - unitCircle)/2;
    unitCircleX = [ unitCircleX; cos(unitCircle) ];
    unitCircleY = [ unitCircleY; sin(unitCircle) ];
    transformedUnitCircleX = [ transformedUnitCircleX;   (1 - (2*cos(psi)).^b.*cos(b*psi)) ];
    transformedUnitCircleY = [ transformedUnitCircleY;   (    (2*cos(psi)).^b.*sin(b*psi)) ];
    
    % Plot the unit circle and its image under the mapping
    % along with the roots of the characterisitc polynomial.
    plot(transformedUnitCircleX, transformedUnitCircleY);
    hold on;
    plot(unitCircleX, unitCircleY);
    scatter(real(charPolyRoots), imag(charPolyRoots));
    maxXYaxis = max( [ transformedUnitCircleX; unitCircleX;...
                       transformedUnitCircleY; unitCircleY ] );
    minXYaxis = min( [ transformedUnitCircleX; unitCircleX;...
                       transformedUnitCircleY; unitCircleY ] );
    maxXYaxis = max( maxXYaxis, -minXYaxis );
    axis( 2*[ -maxXYaxis; maxXYaxis; -maxXYaxis; maxXYaxis; ] );
    hold off;
    title('Roots of the Characteristic Polynomial with the Image of the Unit Circle');
    pause
    
end


% end