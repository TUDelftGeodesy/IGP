function  [B,d]=symdiags(A,d)
%SYMDIAGS  Diagonals from symmetric matrix.
%   [B,d]= SYMDIAGS(A) extracts all cyclic diagonals from the m-by-m symmetric 
%   matrix A.  B is a m-by-p matrix whose columns are the p cyclic diagonals
%   of A, with p = m/2+1. The columns of B include the cyclic diagonals
%   d=0:m/2. The first column of B contains the main diagonal (d(1)=0), the 
%   other columns contains the elements of the d(k)'th and m-d(k)'th diagonals. 
%   The first m-d(k) elements of a column in B are from diagonal d(k), the last
%   d(k) elements are from diagonal d(m-k).  
%
%   B= SYMDIAGS(A,d) extracts the p cyclic diagonals specified by d. The 
%   values of d must be in the range [0:m/2]
%
%   A= SYMDIAGS(B) reconstucts the symmetric m-by-m matrix A from the 
%   m-by-m/2+1 matrix of cyclic diagonals B.  The output matrix is
%   a full matrix.
%
%   A= SYMDIAGS(B,d) reconstucts the symmetric m-by-m matrix A from the p 
%   cyclic diagonals in the m-by-p matrix B, with the cyclic diagonals 
%   specified in d. The values of d must be in the range [0:m/2].
%   The output is sparse matrix with only the upper triangle with
%   diagonals d(1:p) and m-d(1:p).
%
%   Roughly, A, B and d are related by
%
%     for k = 1:p
%        B(:,k) = [ diag(A,d(k)) ; diag(A,m-d(k)) ]
%     end
%
%   See also: DIAG, SPDIAGS
%
%   (c) Hans van der Marel, Delft University of Technology, 2020

%   Created:  11 August 2020 by Hans van der Marel
%   Modified:

% Check the number of input arguments

if nargin < 1 || nargin > 2
    error('incorrect number of input arguments')
end

% Check the input arguments and provide default values

[m,n] = size(A);

kmax=floor(m/2);
if nargin < 2
   d=0:kmax;
end
if any(d < 0)
    error('diagonal d must be positive numbers') 
end
if any(d > kmax )
    error('diagonal d must be smaller than m/2') 
end
p=length(d);  


if m == n

   % Input matrix A is square --> form matrix of cyclic  diagonals B.

   B = zeros(m,p,class(A));
   for k = 1:p
      % store diagonal d(k) in first m-d(k) elements
      B(1:m-d(k),k) = diag(A,d(k));
      % store diagonal m-d(k) in last d(k) elements
      if d(k) > 0 && ~( d(k) == m/2 && d(k)*2 == m )
         B(m-d(k)+1:m,k) = diag(A,m-d(k));
      end
      % please note that when m is even, the last diagonal would be stored
      % twice, which we don't want and prevent with the test d(k)==m/2 && d(k)*2==m
   end
   
else

   % Input matrix A is matrix of cyclic diagonals --> form square matrix

   % Check size of matrix A (should be m-by-p)
   if n ~= p
       error('Incorrect number of columns in input matrix A with cyclic diagonals') 
   end
   
   % The output is a full matrix if d is not specified, if d is present,
   % output is a sparse matrix, with only the upper triangle
   
   sparseOut=true;
   if nargin < 2
      sparseOut=false;
   end
   
   % Compute indices and store values in a compact format row-index,column-index,value
   
   if sparseOut
      b=zeros(m*n,3);
   else
      B = zeros(m,m,class(A));   
   end
   
   for k = 1:p
      % Process the d(k)'th and m-d(k)'th  diagonals
      i=1:m-d(k);
      j(i)=i+d(k);
      if d(k) > 0 
         i(m-d(k)+1:m)=1:d(k);
         j(m-d(k)+1:m)=i(m-d(k)+1:m)+m-d(k);
      end
      if sparseOut
         b((k-1)*m+1:k*m,:) = [i(:) j(:) A(:,k)];
      else       
         if d(k) == m/2 && d(k)*2 == m 
            % please note that when m is even we have a special case for the
            % last diagonal, because of zeros in the input matrix 
            B(sub2ind([m m],i(1:d(k)),j(1:d(k))))=A(1:d(k),k);
         else
            B(sub2ind([m m],i(:),j(:)))=A(:,k);
         end
      end
   end
   
   if sparseOut
      % make sparse output matrix with only the upper triangle
      B = sparse(b(:,1),b(:,2),b(:,3),m,m);
   else
      % make full symmetric matrix
      B=B+triu(B,1)';
   end
      
end

    
end