





  
  # ----- CHECK RANK CONDITION ---- %
    
    p1 <- p + opt$rConstant
  rankJ <- [] %initialize the rank
  
  if(r>0) # If rank is zero then Alpha and Beta are empty
  
  if(is.null(opt$R_Beta))
    H_beta <- eye(p1*r)
  else
    H_beta <- null(opt$R_Beta)
  end
  
  # We use the commutation matrix K_pr to transform vec(A) into vec(A'), # '
  #   see Magnus & Neudecker (1988, p. 47, eqn (1)).
  Ip <- eye(p)
  Kpr <- reshape(kron(Ip(:),eye(r)),p*r,p*r)
  if(is.null(opt$R_Alpha))
  A <- Kpr * eye(p*r)
  else
  A <- null(opt$R_Alpha*inv(Kpr))
  end
  
  rA <- size(A,2) # number of free parameters in alpha
  rH <- size(H_beta,2) # number of free parameters in beta (including constant)
  
  # Following Boswijk & Doornik (2004, p.447) identification condition
  kronA <- kron(eye(p), [results$coeffs.betaHat results$coeffs.rhoHat])...
  *[A zeros(p*r,rH)]
  kronH <- kron(results$coeffs.alphaHat, eye(p1))*[zeros(p1*r,rA) H_beta]
  rankJ <- rank(kronA + kronH)
  
  results$rankJ <- rankJ
  
  end
  
  
  # --- CHECK RANKS OF ALPHA AND BETA ---%
  
  # Check that alpha and beta have full rank to ensure that restrictions
  #   do not reduce their rank.
  if(rank(results$coeffs.alphaHat) < r)
  fprintf('\nWarning: Alpha hat has rank less than r!\n')
  end
  
  if( rank(results$coeffs.betaHat) < r)
  fprintf('\nWarning: Beta hat has rank less than r!\n')
  end
  
  
  # --- FREE PARAMETERS --- %
  
  # Compute the number of free parameters in addition to those in alpha
  #   and beta.
  [ fp ] <- FreeParams(k, r, p, opt, rankJ)
  # Store the result.
  results$fp <- fp
  
  
  # --- STANDARD ERRORS --- %
  
  if(opt$CalcSE)
  # If any restrictions have been imposed, the Hessian matrix must be
  #   adjusted to account for them.
  if( !is.null(opt$R_Alpha) ||!is.null(opt$R_psi) )
  
  # Create R matrix with all restrictions.
  
  # Count the number of restrictions on d,b. Note: opt$R_psi already
  #  contains restrict DB, so the size() is only reliable if
  #  it's turned off. %'
  if(!is.null(opt$R_psi))
  if(opt$restrictDB)
  rowDB <- size(opt$R_psi,1) - 1
  else
  rowDB <- size(opt$R_psi,1)
  end
  else
  # Otherwise d,b are unrestricted.
  rowDB <- 0
  end
  
  # Number of restrictions on alpha.
  rowA  <- size(opt$R_Alpha,1)
  
  # Count the variables.
  colDB <- 1 + !opt$restrictDB
  colA <- p*r
  colG <- p*p*k
  colMu <- opt$levelParam*p
  colRh <- opt$unrConstant*p
  # Length of vec(estimated coefficients in Hessian).
  R_cols  <- colDB + colMu + colRh + colA + colG
  
  # The restriction matrix will have rows equal to the number of
  #   restrictions.
  R_rows <- rowDB + rowA
  
  R <- zeros(R_rows, R_cols)
  
  # Fill in the matrix R.
  
  # Start with restrictions on (d,b) and note that if the model
  #   with d=b is being estimated, only the first column of R_psi is
  #   considered. In that case, if there are zeros found in the first
  #   column, the user is asked to rewrite the restriction.
  if(rowDB>0)
  if(opt$restrictDB)
  # If the model d=b is being estimated, only one
  #  restriction can be imposed and that is on d.
  R(1:rowDB,1:colDB) <- opt$R_psi(1,1)
  else
  R(1:rowDB,1:colDB) <- opt$R_psi
  end
  end
  # Put the R_Alpha matrix into the appropriate place in R.
  if(!is.null(opt$R_Alpha))
  R(1+rowDB: rowDB + rowA, ...
  1 + colDB + colMu + colRh: colDB + colMu + colRh + colA)...
  <- opt$R_Alpha
  end
  
  # Calculate unrestricted Hessian.
  [ H ] <- FCVARhess(x, k, r, results$coeffs, opt)
  
  
  # Calculate the restricted Hessian.
  Q <- -inv(H) + inv(H)*R'*inv(R*inv(H)*R')*R*inv(H)
                              else
                              # Model is unrestricted.
                              [ H ] <- FCVARhess(x, k, r, results$coeffs, opt)
                              Q <- -inv(H)
                              end
                              else
                              NumCoeffs <- length(SEmat2vecU(results$coeffs, k, r, p, opt))
                              Q <- zeros(NumCoeffs)
                              end
                              # Calculate the standard errors and store them.
                              SE <- sqrt(diag(Q))
                              results$SE <- SEvec2matU(SE,k,r,p, opt)
                              results$NegInvHessian <- Q
                              
                              
                              
                              # --- GET RESIDUALS --- %
                              
                              [ epsilon ] <- GetResiduals(x, k, r, results$coeffs, opt)
                              results$Residuals <- epsilon
                              
                              
                              # --- OBTAIN ROOTS OF CHARACTERISTIC POLYNOMIAL --- %
                              
                              cPolyRoots <- CharPolyRoots(results$coeffs, opt, k, r, p)
                              results$cPolyRoots <- cPolyRoots
                              
                              
                              # --- PRINT OUTPUT --- %
                              
                              if (opt$print2screen)
                              if(!opt$CalcSE)
                              fprintf('Warning: standard errors have not been calculated!\n')
                              end
                              
                              # create a variable for output strings
                              yesNo <- {'No','Yes'}
                              fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n')
                              fprintf(1,'                      Fractionally Cointegrated VAR: Estimation Results                              ')
                              fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n')
                              fprintf(1,'Dimension of system:  %6.0f      Number of observations in sample:       %6.0f \n', p, T+opt$N)
                              fprintf(1,'Number of lags:       %6.0f      Number of observations for estimation:  %6.0f \n', k, T)
                              fprintf(1,'Restricted constant:  %6s      Initial values:                         %6.0f\n', yesNo{opt$rConstant+1}, opt$N )
                              fprintf(1,'Unrestricted constant:%6s      Level parameter:                        %6s\n', yesNo{opt$unrConstant+1}, yesNo{opt$levelParam+1} )
                              if(size(opt$R_psi,1)==1)
                              # 1 restriction.
                              dbUB <- H_psi*UB(1)
                              dbLB <- H_psi*LB(1)
                              dbStart <- H_psi*startVals(1)
                              fprintf(1,'Starting value for d:    %1.3f    Parameter space for d: (%1.3f , %1.3f) \n', dbStart(1), dbLB(1), dbUB(1))
  fprintf(1,'Starting value for b:    %1.3f    Parameter space for b: (%1.3f , %1.3f) \n', dbStart(2), dbLB(2), dbUB(2))
else
# Unrestricted or 2 restrictions.
fprintf(1,'Starting value for d:    %1.3f    Parameter space for d: (%1.3f , %1.3f) \n', startVals(1), LB(1), UB(1))
fprintf(1,'Starting value for b:    %1.3f    Parameter space for b: (%1.3f , %1.3f) \n', startVals(2), LB(2), UB(2))
fprintf(1,'Imposing d >= b:      %6s\n', yesNo{opt$constrained+1} )

end
fprintf(1,'-----------------------------------------------------------------------------------------------------\n')
fprintf(1,'Cointegrating rank:   %10.0f  AIC:            %10.3f \n', r, -2*maxLike + 2*fp)
fprintf(1,'Log-likelihood:       %10.3f  BIC:            %10.3f \n', maxLike, -2*maxLike + fp*log(T))
fprintf(1,'log(det(Omega_hat)):  %10.3f  Free parameters:%10.0f \n', log(det(results$coeffs.OmegaHat)), fp)
fprintf(1,'-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '    Fractional parameters:                                                                             \n')
fprintf(1,    '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '    Coefficient              \t Estimate              \t  Standard error \n')
fprintf(1,    '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '         d                   \t %8.3f              \t     %8.3f                \n', results$coeffs.db(1), results$SE.db(1))
if !opt$restrictDB
fprintf(1,'         b                   \t %8.3f              \t     %8.3f                \n', results$coeffs.db(2), results$SE.db(2))
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
if r > 0
if opt$rConstant
varList <- '(beta and rho):'
else
varList <- '(beta):        '
end
fprintf(1,'    Cointegrating equations %s                                                          \n', varList)
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '      Variable      ' )
for j <- 1:r
fprintf(1,    '  CI equation %d  ', j)
end
fprintf(1,'\n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for i <- 1:p
fprintf(1,    '        Var%d       ',i )
for j <- 1:r
fprintf(1,'    %8.3f     ', results$coeffs.betaHat(i,j) )
end
fprintf(1,'\n')
end
if opt$rConstant
fprintf(1,    '      Constant     ' )
for j <- 1:r
fprintf(1,'    %8.3f     ', results$coeffs.rhoHat(j) )
end
fprintf(1,'\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
if (is.null(opt$R_Alpha) && is.null(opt$R_Beta) )
fprintf(1,  'Note: Identifying restriction imposed.                                                               \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
end
fprintf(1,'    Adjustment matrix (alpha):                                                                         \n' )
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '      Variable      ' )
for j <- 1:r
fprintf(1,    '  CI equation %d  ', j)
end
fprintf(1,'\n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for i <- 1:p
fprintf(1,    '        Var %d      ',i )
for j <- 1:r
fprintf(1,'    %8.3f     ', results$coeffs.alphaHat(i,j) )
end
fprintf(1,'\n')
fprintf(1,    '         SE %d      ',i )
for j <- 1:r
fprintf(1,'   (%8.3f  )  ', results$SE.alphaHat(i,j) )
end
fprintf(1,'\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  'Note: Standard errors in parenthesis.                                                                \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,'    Long-run matrix (Pi):                                                                       \n' )
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '      Variable  ' )
for j <- 1:p
fprintf(1,    '       Var %d   ', j)
end
fprintf(1,'\n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for i <- 1:p
fprintf(1,    '      Var %d      ',i )
for j <- 1:p
fprintf(1,'   %8.3f    ', results$coeffs.PiHat(i,j) )
end
fprintf(1,'\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
end
end

# Print level parameter if present.
if (opt$print2screen && opt$levelParam)
fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n')
fprintf(1,'    Level parameter (mu):                                                                         \n' )
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for i <- 1:p
fprintf(1,    '        Var %d      ',i )
fprintf(1,'    %8.3f     ', results$coeffs.muHat(i) )
fprintf(1,'\n')
fprintf(1,    '         SE %d      ',i )
fprintf(1,'   (%8.3f  )  ', results$SE.muHat(i) )
fprintf(1,'\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  'Note: Standard errors in parenthesis (from numerical Hessian) but asymptotic distribution is unknown. \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
end

# Print unrestricted constant if present.
if (opt$print2screen && opt$unrConstant)
fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n')
fprintf(1,'    Unrestricted constant term:                                                                     \n' )
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for i <- 1:p
fprintf(1,    '        Var %d      ',i )
fprintf(1,'    %8.3f     ', results$coeffs.xiHat(i) )
fprintf(1,'\n')
fprintf(1,    '         SE %d      ',i )
fprintf(1,'   (%8.3f  )  ', results$SE.xiHat(i) )
fprintf(1,'\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  'Note: Standard errors in parenthesis (from numerical Hessian) but asymptotic distribution is unknown. \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
end

# Print Gamma coefficients if required.
if opt$print2screen && opt$printGammas && (k > 0)
for l <- 1:k
GammaHatk <- results$coeffs.GammaHat( :, p*(l-1)+1 : p*l )
GammaSEk <- results$SE.GammaHat( :, p*(l-1)+1 : p*l )

fprintf(1,'    Lag matrix %d (Gamma_%d):                                                                            \n', l, l )
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,    '      Variable  ' )
for j <- 1:p
fprintf(1,    '       Var %d   ', j)
end
fprintf(1,'\n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for i <- 1:p
fprintf(1,    '      Var %d      ',i )
for j <- 1:p
fprintf(1,'   %8.3f    ', GammaHatk(i,j) )
end
fprintf(1,'\n')
fprintf(1,    '       SE %d       ',i )
for j <- 1:p
fprintf(1,' (%8.3f  )  ', GammaSEk(i,j) )
end
fprintf(1,'\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  'Note: Standard errors in parentheses.                                                                \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
end
end

# Print roots of characteristic polynomial if required.
if (opt$print2screen && opt$printRoots)
fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  '    Roots of the characteristic polynomial                                                           \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  '    Number     Real part    Imaginary part       Modulus                                             \n')
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n')
for j <- 1:length(cPolyRoots)
fprintf(1, '      %2.0f       %8.3f       %8.3f         %8.3f                                        \n',...
j, real(cPolyRoots(j)), imag(cPolyRoots(j)), abs(cPolyRoots(j)) )
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n')
end

# Print notifications regarding restrictions.
if (opt$print2screen &&  (!is.null(opt$R_Alpha) || !is.null(opt$R_psi) ...
|| !is.null(opt$R_Beta) ))
fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n')
fprintf(1,  'Restrictions imposed on the following parameters:\n')
if(!is.null(opt$R_psi))
fprintf('- Psi. For details see "options$R_psi"\n')
end
if(!is.null(opt$R_Alpha))
fprintf('- Alpha. For details see "options$R_Alpha"\n')
end
if(!is.null(opt$R_Beta))
fprintf('- Beta. For details see "options$R_Beta"\n')
end
fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n')
end

end
