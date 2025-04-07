function y=stmtransapply(y,xref,xoff)
%stmtransapply  Apply Space Time Matrix S-Transformation
%  Y=STMTRANSAPPLY(Y,XREF,XOFF) applies the S-Transformation to the 
%  space-time matrix Y with XREF the transforming reference time series 
%  and XOFF the transforming offsets. Y can be a 2D or 3D matrix with
%  dimensions [np,ne,no], XREF and XOFF are a vector or matrix with 
%  dimensions of respectively [no,ne] and [np,no]. 
%
%  The function STMTRANS can be used to compute XREF and XOFF using a
%  variety of methods.
%
%  See also stmtrans.
%
%  (c) Hans van der Marel, Delft University of Technology, 2024.

% Created:   6 May 2024 by Hans van der Marel
% Modified: 

% check dimensions and defaults

[np,ne,no]=size(y);

if isempty(xref)
    xref=zeros(no,ne);
end
if isempty(xoff)
    xoff=zeros(np,no);
end

if size(xref,1)~=no || size(xref,2)~=ne 
   % try tranpose
   xref=xref';
end
assert(size(xref,1)==no && size(xref,2)==ne  ,'stmtransapply dimensions of xref must match y.')
assert(size(xoff,1)==np && size(xoff,2)==no  ,'stmtransapply dimensions of xoff must match y.')

% actual S-transform

for k=1:no
   y(:,:,k) = y(:,:,k) - repmat(xref(k,:),[np,1]) - repmat(xoff(:,k),[1,ne]);
end

end
