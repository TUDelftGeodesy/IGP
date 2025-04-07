function [H,qmax] = mht_setup_models(hyp,Bt,BT,minEpoch)

if size(Bt,1)==1
    Bt = Bt'; %Force column vector
end
if size(BT,1)==1
    BT = BT'; %Force column vector
end

qmax = 0;
Nhyp = length(hyp);
Nepoch = length(Bt);

count = 0;
for v = 1:Nhyp

  models = hyp{v};
  Nmodels = length(models);
  m = [];

  for w = 1:Nmodels

    switch char(models(w))
     case 'constant'
       m(w).model = ones(size(Bt));
       m(w).N = 1;
       m(w).Npar = 1;
       m(w).epoch = NaN;
     case 'linear'
       m(w).model = Bt;
       m(w).N = 1;
       m(w).Npar = 1;
       m(w).epoch = NaN;
     case 'temperature'
       m(w).model = BT;
       m(w).N = 1;
       m(w).Npar = 1;
       m(w).epoch = NaN;
     case 'periodic'
       m(w).model = [sin(2*pi*Bt) cos(2*pi*Bt)-1];
       m(w).N = 1;
       m(w).Npar = 2;
       m(w).epoch = NaN;
     case 'breakpoint'
       temp = repmat(Bt,1,Nepoch);
       temp = temp-temp';
       temp = tril(temp);
       m(w).model = temp(:,minEpoch+1:end-minEpoch);
       m(w).N = Nepoch - 2*minEpoch;
       m(w).Npar = 1;
       m(w).epoch = [minEpoch+1:Nepoch-minEpoch]';
     case 'heaviside'
       temp = tril(ones(Nepoch));
       m(w).model = temp(:,minEpoch+1:end-minEpoch);
       m(w).N = Nepoch - 2*minEpoch;
       m(w).Npar = 1;
       m(w).epoch = [minEpoch+1:Nepoch-minEpoch]';
     otherwise
       error('You specified an unknown model');
     end

  end

  if m(1).Npar>1 & m(1).N==1 %
    A = m(1).model;
  else
    A = m(1).model(:); % if m(1).N > 1, matrix is converted to long vector
  end
  epoch = m(1).epoch;
  Nalt = m(1).N;

  for w = 2:Nmodels
    Anew = repmat(m(w).model,Nalt,1);
    if m(w).Npar>1
      A = [kron(ones(m(w).N,1),A) Anew];
    else
      A = [kron(ones(m(w).N,1),A) Anew(:)];
    end
    epoch = [kron(ones(m(w).N,1),epoch) kron(m(w).epoch,ones(Nalt,1))];
    Nalt = Nalt*m(w).N;
  end

  q = size(A,2); %dimension of model
  qmax = max(q,qmax);
  for w = 1:Nalt
    count = count+1;
    H(count).models = models;
    H(count).q = q;
    H(count).Cj = A((w-1)*Nepoch+1:w*Nepoch,:);
    H(count).epoch = epoch(w,:);
  end

  clear q A Anew Nalt epoch

end

