%% Long Random (de-)interleaver

%% Inputs:
%  BitsIn   -sequence of 1 x Nb bit sequence to be interleaved
%  L        -interleaving length (in number of bits)      
%  varargin{1} - matrix round(Nb/L) x L representing interleaving maps
%  (when working as a deinterleaver)

%% Outputs:
%  varargout{1} - (de)/interleaved bit sequence
%  varargout{2} - matrix round(Nb/L) x L representing interleaving maps.
%  This is randomly generated when the function is used as an interleaver

function varargout = LongRandInter(BitsIn,L,varargin)
if nargin==2 % Works as an interleaver
varargout{1}=BitsIn;
r=mod(length(BitsIn),L);    
q=(length(BitsIn)-r)/L;
map=cell(1,q);
for qq=1:q
if qq<q
map{qq}=randperm(L);
idx=1+(qq-1)*L:qq*L;
else
map{qq}=randperm(L+r);
idx=1+(qq-1)*L:qq*L+r;   
end
Chop=BitsIn(idx);
varargout{1}(idx)=Chop(map{qq});
end
varargout{2}=map;

elseif nargin>2 % Works as a deinterleaver
map=varargin{1};
BitsOut=BitsIn;
q=length(map);
L=length(map{1});
for qq=1:q
  if qq<q
  idx=1+(qq-1)*L:qq*L;
  else
  r=mod(length(BitsIn),L); 
  idx=1+(qq-1)*L:qq*L+r;   
  end
if nargin==3   
Chop=BitsIn(idx);
BitsOut(map{qq}+(qq-1)*L)=Chop; 
elseif nargin==4 && varargin{2}==-1
BitsOut(idx)=BitsIn(map{qq}+(qq-1)*L);    
end
varargout{1}=BitsOut;
end

else 
  error('Wrong number of input arguments');
end

end