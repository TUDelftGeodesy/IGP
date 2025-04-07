function qout=covreformat2(qin,fmtin,fmtout) 
%COVREFORMAT  Reformat co-variance matrix.
%   QOUT=COVREFORMAT(QIN,FMTIN,FMTOUT) reformats the co-variance matrix
%   QIN with format FMTIN into QOUT with format FMTOUT. FMTIN and FMTOUT
%   are strings with the format identifier. 
%
%   Supported formats FMTIN and FMTOUT are:
%
%     qmat    co-variance matrix [ qx qxy qxz ; qxy qy qyz ; qxy qyz qz ]
%
%     qvec    vector with qx, qy, qz, qxy, qxz, qyz 
%     scor    vector with sx, sy, sz, rxy, rxz, ryz     (geo++, NRCAN)
%     scof    vector with sx, sy, sz, cxy, cxz, cyz     (gamit/globk)
%     scov    vector with sx, sy, sz, sxy, sxz, syz
%
%     qvecd   vector with qx, qy, qz, qxy, qyz, qxz 
%     scord   vector with sx, sy, sz, rxy, ryz, rxz     
%     scofd   vector with sx, sy, sz, cxy, cyz, cxz     
%     scovd   vector with sx, sy, sz, sxy, syz, sxz     (rtklib)
%
%   with si=sqrt(qi), rij=qij/(sqrt(qi)*sqrt(qj)),
%   cij=sign(rij)*sqrt(abs(rij)) and sij=sign(qij)*sqrt(abs(qij)).
%
%   QIN and QOUT can have the following shapes: Q(1:3,1:3) or q(1:6), or, 
%   Q(1:3,1:3,k) or q(k,1:6) for k=1...N.
%
%   This function generates a warning in case the co-variance matrix is
%   not positive definite and substitutes the matrix with NaN's.
%
%   Example:
%
%      scov = [ 0.2061  0.0708  0.6352  -0.0109  0.1514  -0.0188 ]
%      Q=covreformat(scov,'scovd','qmat')
%
%   (c) Hans van der Marel, Delft University of Technology, 2014

%   Created:     8 March 2014 by Hans van der Marel
%   Modified:   11 July 2018 by Hans van der Marel
%                  - new formats scof and scofd (square root of correlation)

% Check input arguments

squeze=false;
if strcmpi(fmtin,'qmat') 
   if size(qin,1) ~= 3 || size(qin,2) ~= 3
      error('input covariance matrix must be 3-by-3')
   end
   if ndims(qin) == 2
     squeze=true;
     qin(:,:,1)=qin(:,:);
     n=1;
   else
     n=size(qin,3);
   end
else
   if ndims(qin) == 1
     squeze=true;
     qin(1,:)=qin;
     n=1;
   else
     n=size(qin,1);
   end
   if size(qin,2) ~= 6 
      error('input covariance matrix must be a vector or matrix with 6 columns')
   end     
end

% Preallocate memory

switch lower(fmtout)
    case 'qmat'
       qout=zeros(3,3,n);
    otherwise
       qout=zeros(n,6);
end

% Convert scof and scofd format to scor or scord

if strncmpi(fmtin,'scof',4)
   qin(:,4:6)=sign(qin(:,4:6)).*abs(qin(:,4:6)).^2;
end

% For each set convert

for k=1:n
    
  % First convert input format to plain co-variance matrix qxyz

  switch lower(fmtin)
    case 'qmat'
        qxyz=qin(:,:,k);
    case 'qvec'
        qxyz= [ qin(k,1)  qin(k,4)  qin(k,5) ; ...
                qin(k,4)  qin(k,2)  qin(k,6) ; ...
                qin(k,5)  qin(k,6)  qin(k,3) ];
    case 'qvecd'
        qxyz= [ qin(k,1)  qin(k,4)  qin(k,6) ; ...
                qin(k,4)  qin(k,2)  qin(k,5) ; ...
                qin(k,6)  qin(k,5)  qin(k,3) ];
    case {'scor','scof'}
        % sdx sdy sdz rxy rxz ryz
        qxyz= [ qin(k,1)^2                   qin(k,4)*(qin(k,1)*qin(k,2))   qin(k,5)*(qin(k,1)*qin(k,3)) ; ...
                qin(k,4)*(qin(k,2)*qin(k,1)) qin(k,2)^2                     qin(k,6)*(qin(k,2)*qin(k,3)) ; ...
                qin(k,5)*(qin(k,3)*qin(k,1)) qin(k,6)*(qin(k,3)*qin(k,2))   qin(k,3)^2 ];
    case {'scord','scofd'}
        % sdx sdy sdz rxy ryz rxz
        qxyz= [ qin(k,1)^2                   qin(k,4)*(qin(k,1)*qin(k,2))   qin(k,6)*(qin(k,1)*qin(k,3)) ; ...
                qin(k,4)*(qin(k,2)*qin(k,1)) qin(k,2)^2                     qin(k,5)*(qin(k,2)*qin(k,3)) ; ...
                qin(k,6)*(qin(k,3)*qin(k,1)) qin(k,5)*(qin(k,3)*qin(k,2))   qin(k,3)^2 ];
    case 'scov'
        % sdx sdy sdz sdxy sdxz sdyz
        qxyz= [ qin(k,1)^2                 sign(qin(k,4))*qin(k,4)^2  sign(qin(k,5))*qin(k,5)^2 ; ...
                sign(qin(k,4))*qin(k,4)^2  qin(k,2)^2                 sign(qin(k,6))*qin(k,6)^2 ; ...
                sign(qin(k,5))*qin(k,5)^2  sign(qin(k,6))*qin(k,6)^2  qin(k,3)^2                  ];
    case 'scovd'
        % sdx sdy sdz sdxy sdyz sdxz
        qxyz= [ qin(k,1)^2                 sign(qin(k,4))*qin(k,4)^2  sign(qin(k,6))*qin(k,6)^2 ; ...
                sign(qin(k,4))*qin(k,4)^2  qin(k,2)^2                 sign(qin(k,5))*qin(k,5)^2 ; ...
                sign(qin(k,6))*qin(k,6)^2  sign(qin(k,5))*qin(k,5)^2  qin(k,3)^2                  ];
     otherwise
        error('unknown input covariance matrix format')
  end
  
  % Check that the co-variance matrix is positive definite
  
  try
    chol(qxyz); 
  catch exception
    disp('co-variance matrix is not positive definite')
    qxyz=nan(3,3);
  end
  
  % Next convert plain co-variance matrix qxyz into output format

  switch lower(fmtout)
    case 'qmat'
        qout(:,:,k)=qxyz(:,:);
    case 'qvec'
        qout(k,:)=[ qxyz(1,1) qxyz(2,2) qxyz(3,3) qxyz(1,2) qxyz(1,3) qxyz(2,3) ];
    case 'qvecd'
        qout(k,:)=[ qxyz(1,1) qxyz(2,2) qxyz(3,3) qxyz(1,2) qxyz(2,3) qxyz(1,3) ];
    case {'scor','scof'}
        % sdx sdy sdz rxy rxz ryz
        qout(k,:)=[ sqrt(qxyz(1,1)) sqrt(qxyz(2,2)) sqrt(qxyz(3,3)) ...
                    qxyz(1,2)/sqrt(qxyz(1,1)*qxyz(2,2))  qxyz(1,3)/sqrt(qxyz(1,1)*qxyz(3,3)) qxyz(2,3)/sqrt(qxyz(2,2)*qxyz(3,3)) ];
    case {'scord','scofd'}
        % sdx sdy sdz rxy ryz rxz
        qout(k,:)=[ sqrt(qxyz(1,1)) sqrt(qxyz(2,2)) sqrt(qxyz(3,3)) ...
                    qxyz(1,2)/sqrt(qxyz(1,1)*qxyz(2,2))  qxyz(2,3)/sqrt(qxyz(2,2)*qxyz(3,3)) qxyz(1,3)/sqrt(qxyz(1,1)*qxyz(3,3)) ];
    case 'scov'
        % sdx sdy sdz sdxy sdxz sdyz
        qout(k,:)=[ sqrt(qxyz(1,1)) sqrt(qxyz(2,2)) sqrt(qxyz(3,3)) ...
               sign(qxyz(1,2))*sqrt(abs(qxyz(1,2)))  sign(qxyz(1,3))*sqrt(abs(qxyz(1,3))) sign(qxyz(2,3))*sqrt(abs(qxyz(2,3))) ];
    case 'scovd'
        % sdx sdy sdz sdxy sdyz sdxz
        qout(k,:)=[ sqrt(qxyz(1,1)) sqrt(qxyz(2,2)) sqrt(qxyz(3,3)) ...
               sign(qxyz(1,2))*sqrt(abs(qxyz(1,2)))  sign(qxyz(2,3))*sqrt(abs(qxyz(2,3))) sign(qxyz(1,3))*sqrt(abs(qxyz(1,3))) ];
     otherwise
        error('unknown output covariance matrix format')
  end

end

% Optionally convert scor or scord output format to scof or scofd

if strncmpi(fmtout,'scof',4)
   qout(:,4:6)=sign(qout(:,4:6)).*sqrt(abs(qout(:,4:6)));
end

% Reproduce the same output format

if squeze && strcmpi(fmtout,'qmat') 
   qout(:,:)=qout(:,:,k);
elseif ndims(qin) == 1
   qout=qout(1,:);
end

end