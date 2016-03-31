# Adopted from Matt Landreman's Chebyshev construction routines in Matlab (2014)
module chebyshev

export singleChebyshevWeights, singleChebyshevDifferentiation

function singleChebyshevWeights(N1, a, b)
  N=N1-1
  bma=b-a
  c=zeros(N1,2)
  c[1:2:N1,1]=(2./[1; 1-(2:2:N).^2 ])'
  c[2,2]=1

  f=[real(ifft([c[1:N1,1]; c[N:-1:2,1]])) real(ifft([c[1:N1,2];c[N:-1:2,2]]))]
  w=bma*([f[1,1]; 2*f[2:N,1]; f[N1,1]])/2
  x= 0.5 *((b+a)+N*bma*f[1:N1,2])
  x=x[end:-1:1]
  w=w[end:-1:1]'
  return x,w
end

function singleChebyshevDifferentiation(N, xMin, xMax)
  N1=N-1
  if N1==0
     D=0
     x=1
     return x,D
  end
  x = cos(pi*[0:N1]/N1)
  c = [2; ones(N1-1,1); 2].*(-1).^(0:N1)
  X = repmat(x,1,N1+1)
  dX = X-X'
  D = (c*(1./c)')./(dX+(eye(N1+1)))
  
  D = D - diagm(vec(sum(D',1)))
  D = D * 2/(xMax-xMin);
  x = (x+1) * (xMax-xMin)/2 + xMin
  D=fliplr(flipud(D));
  x=fliplr(x');
  return x,D
end


end
