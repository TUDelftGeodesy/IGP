function [yq,syq,nyq]=mainterp1(x,y,xq,xwindow)
%MMINTERP1  Moving Average interpolation
%   [YQ,SYQ,NYQ]=MAINTERP1(X,Y,XQ) computes the moving average YQ at data
%   points XQ from the input data series Y with abcissae X. SYQ is the 
%   the standard deviation. NYQ is an array with the number of
%   data points that have been used. XQ must be increasing.
%
%   [...]=MMFILT(...,XWINDOW) uses a fixed window length of XWINDOW.
%   Default is to define borders at (XQ(K)+XQ(K+1))/2. 
%
%   (c) TU Delft, Hans van der Marel, 2018

if nargin < 3, error('Need to specify at least X, Y and XQ'); end
if nargin > 4, error('Too many input arguments, please check'); end

if size(x,1) == 1
  x=x';
  y=y';
end
xq=xq(:);

n=size(y,1);
if size(x,1)~=n, 
    error('X and Y must have the same number of columns'); 
end

nq=numel(xq);
if nargin < 4
   xq2=median(diff(xq(:)))/2;
   xqq(1,1)=xq(1)-xq2;
   xqq(2:nq,1)=(xq(1:nq-1)+xq(2:nq))/2;
   xqq(1:nq-1,2)=xqq(2:end,1);
   xqq(nq,2)=xq(nq)+xq2;
else
   xq2=xwindow/2;
   xqq(:,1)=xq-xq2;
   xqq(:,2)=xq+xq2;
end

yq=nan(size(xq,1),size(y,2));
syq=nan(size(xq,1),size(y,2));
nyq=nan(size(xq));

for k=1:nq
   idx=x > xqq(k,1) & x <= xqq(k,2);
   yq(k,:)=nanmean(y(idx,:),1);
   nyq(k)=sum(idx);
   syq(k,:)=nanstd(y(idx,:),0,1);
end

end