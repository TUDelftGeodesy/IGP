function [xref,xoff,y]=stmtrans(y,varargin)
%stmtrans      Space Time Matrix S-Transformation.
%  [XREF,XOFF,Y]=STMTRANS(Y) computes a S-Transformation basis for the
%  space-time matrix Y. Y can be a 2D or 3D space time matrix with 
%  dimensions [np,ne,no]. The output is XREF with the transforming reference 
%  time series and XOFF with the transforming offsets, and optionally the
%  S-transformed space-time matrix Y. XREF and XOFF are vectors or matrices 
%  with dimensions [no,ne] and [np,no]. 
%
%  [XREF,XOFF,Y]=STMTRANS(Y, option, value, ...) with options
%
%    method       S-Transform method, possible values
%      'none'                   - do nothing, XREF and XOFF contain zeros
%      'offsets-only','offsets' - only do offsets
%      'minimum-norm','minnorm' - minumum norm solution
%      'sd-norm','sdnorm'       - minimum norm of single differences in 
%                                 time (default)
%      'ssd-norm','ssdnorm'     - minimum norm of symmetric single differences  
%                                 in time
%
%    pntRefMask   Logical array defining subset of points with near zero
%                 displacements (default empty)
%
%  The transform is a one or two step procedure. In the first step the 
%  S-basis is computed using the method given by the 'method' option using
%  all points.
% 
%  Then in the second step, but only if 'pntRefMask' is not empty, corrections 
%  to the reference time series is computed such that the displacements 
%  for points in logical array 'pntRefMask' have minimum-norm. If an 
%  epoch cannot be computed in this way, the correction is interpolated
%  from other epochs.
%
%  The function STMTRANSAPPLY can be used to apply the S-transformation
%  given by XREF and XOFF at a later time. 
%
%  See also stmtransapply.
%
%  (c) Hans van der Marel, Delft University of Technology, 2024.

% Created:   6 May 2024 by Hans van der Marel
% Modified: 16 May 2024 by Hans van der Marel
%              - improved documentation
%              - added some extra checks and warning/error messages
%              - more sensible default options
%           11 June 2024 by Hans van der Marel
%              - build in upper limit to number of epochs for sdnorm methods

% Check the inputs and options

if nargin < 1 
    error('stmtrans expects at least one argument')
end 

opt.method='sd-norm'; 
opt.pntRefMask=[];  
opt.year=[];

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       error(['stmtrans invalid option ' varargin{k} ', abort...'])
    end
end

% Get dimensions and prepare some of the outputs

[np,ne,no]=size(y);
xref=nan(no,ne);
xoff=nan(np,no);

% Check the number of epochs for certain methods

if any(strcmpi(opt.method,{'sd-norm','sdnorm','ssd-norm','ssdnorm'})) && ne > 1000
    warning('stmtrans: too many epochs for sdnorm or ssdnorm method, fall back to none.')
    opt.method='none';
end

% Nothing to do

if strcmpi(opt.method,'none')
   xref=zeros(no,ne);
   xoff=zeros(np,no);
   return
end

% Prepare some inputs

pntRefMask=opt.pntRefMask;
if isempty(pntRefMask)
   pntRefMask=true(np,1);
end

year=opt.year;
if isempty(year)
   year=1:ne;
end

% Compute the S-tranformation (loop over observation types)

for k=1:no

    % Initial S-transform using all available points (part I)
    
    switch lower(opt.method)
      case {'offsets-only','offsets'}
         xref(k,:)=zeros(1,ne);
      case {'minimum-norm','minnorm'}
         [xref(k,:),xoff(:,k)]=refseries1(y(:,:,k));
         xref(k,:)=xref(k,:)-median(xref(k,:));
         y(:,:,k) = y(:,:,k) - repmat(xref(k,:),[np,1]);
      case {'sd-norm','sdnorm'}
         xref(k,:)=refseries2(y(:,:,k),year);
         xref(k,:)=xref(k,:)-median(xref(k,:));
         y(:,:,k) = y(:,:,k) - repmat(xref(k,:),[np,1]);
      case {'ssd-norm','ssdnorm'}
         xref(k,:)=refseries3(y(:,:,k),year);
         xref(k,:)=xref(k,:)-median(xref(k,:));
         y(:,:,k) = y(:,:,k) - repmat(xref(k,:),[np,1]);
      otherwise
         error(['stmtrans invalid method ' opt.method]);
   end
   xoff(:,k)=nanmean(y(:,:,k),2);
   y(:,:,k) = y(:,:,k) - repmat(xoff(:,k),[1,ne]);

   % Optionally refine the transformation using only selected (stable) points

   if ~isempty(opt.pntRefMask) && any(pntRefMask)
      % Update it using a sub-selection of points given in prtRefMask
      pntRefMask1 = pntRefMask & sum(~isnan(y(:,:,k)),2) > 1;
      [xrefupd,xoffupd]=refseries1(y(pntRefMask1,:,k));
      xrefupd=refseries2(y(pntRefMask1,:,k));
      ip=find(isnan(xrefupd));
      if numel(ip) > ne-2
         % less than two epochs, no interpolation possible
         warning('stmtrans less than two epochs with selected reference points, cannot update S-transform with reference points')
         xrefupd=zeros(1,ne);
      elseif ~isempty(ip)
         % we need to interpolate some epochs
         x=1:ne;
         lp=~isnan(xrefupd);
         xrefupd(ip) = interp1(x(lp),xrefupd(lp),ip,'pchip');
      end
      fprintf('stmtrans method %s, with ref-series update using %d (out of %d) points, with %d (out of %d) epochs interpolated.\n',opt.method,sum(pntRefMask1),np,numel(ip),ne)
      disp('number of reference points per epoch:')
      disp(['  ' num2str(sum(~isnan(y(pntRefMask1,:,k)),1))])
      y(:,:,k) = y(:,:,k) - repmat(xrefupd,[np,1]);
      xref(k,:)=xref(k,:)+xrefupd;
      xoff(:,k)=nanmean(y(:,:,k),2);
      y(:,:,k) = y(:,:,k) - repmat(xoff(:,k),[1,ne]);
   else
      fprintf('stmtrans method %s w/o ref-series update.\n',opt.method)
   end

end

end

function [xref,xoff]=refseries1(y)
%refseries1  Minimun-norm reference series solution

[np,ne]=size(y);
m=np*ne;
n=np+ne;

ii=repmat([1:np]',[1,ne]);
jj=repmat([1:ne],[np,1]);
A2=sparse( [1:m 1:m]',[ jj(:); ii(:)+ne],ones(m*2,1),m,n,m*2);

y=double(y(:));
ymask=~isnan(y);

N22=A2(ymask,:)'*A2(ymask,:);

idx=diag(N22) > 0;
%x2=zeros(size(idx));
x2=nan(size(idx));

b2=A2(ymask,:)'*y(ymask);

x2(idx)=N22(idx,idx)\b2(idx);

xref=x2(1:ne)';
xoff=x2(ne+1:end);

end

function xref=refseries2(y,t)
%refseries2  Single difference minimum-norm reference series solution

ne=size(y,2);
if nargin < 2
   t=1:ne;
end

ymat=zeros(ne,ne);
yvec=zeros((ne-1)*ne/2,1);
wvec=zeros((ne-1)*ne/2,1);
Amat=zeros((ne-1)*ne/2,ne);

kk=0;
for k=2:ne
    for l=1:k-1
        ymat(k,l)=nanmean(y(:,k)-y(:,l));
        kk=kk+1;
        yvec(kk)=ymat(k,l);
        wvec(kk)=sum(~isnan(y(:,k)-y(:,l)))/(k-l)^2;
        Amat(kk,k)=1;
        Amat(kk,l)=-1;
    end
end

Amat(isnan(yvec),:)=[];
wvec(isnan(yvec),:)=[];
yvec(isnan(yvec),:)=[];

Nmat=Amat'*diag(wvec)*Amat;
bvec=Amat'*diag(wvec)*yvec;

idx=diag(Nmat) > 0;

xref=nan(size(idx));
xref(idx)=pinv(Nmat(idx,idx))*bvec(idx);
xref=xref';

end

function xref=refseries3(y,t)
%refseries3  Symmetric single difference minimum-norm reference series solution

ne=size(y,2);
if nargin < 2
   t=1:ne;
end

nbin = ne*(ne-1)*(ne-2)/6;

yvec=zeros(nbin,1);
wvec=zeros(nbin,1);
Amat=zeros(nbin,ne);

kk=0;
for k1=1:ne-2
   for k2=k1+1:ne-1
      for k3=k2+1:ne
         kk=kk+1;
         cond = (t(k2)-t(k1))*y(:,k3) + (t(k1)-t(k3))*y(:,k2) + (t(k3)-t(k2))*y(:,k1);
         yvec(kk)=nanmean(cond);
         wvec(kk)=sum(~isnan(cond))/(abs((t(k2)-t(k1))*(t(k1)-t(k3))*(t(k3)-t(k2)))^2);
         %wvec(kk)=sum(~isnan(cond));
         Amat(kk,k3)=t(k2)-t(k1);
         Amat(kk,k2)=t(k1)-t(k3);
         Amat(kk,k1)=t(k3)-t(k2);
      end
   end
end

Amat(isnan(yvec),:)=[];
wvec(isnan(yvec),:)=[];
yvec(isnan(yvec),:)=[];

Nmat=Amat'*diag(wvec)*Amat;
bvec=Amat'*diag(wvec)*yvec;

idx=diag(Nmat) > 0;

xref=nan(size(idx));
xref(idx)=pinv(Nmat(idx,idx))*bvec(idx);
xref=xref';

%diff(xref)
%size(Nmat)
%rank(Nmat)

v=nanmedian(diff(xref)./diff(t));

xref=xref-(t(:)'-t(floor(ne/2)))*v;
xref=xref-nanmedian(xref);

end
