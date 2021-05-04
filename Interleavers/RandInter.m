%% Random (de-)interleaver

function varargout = RandInter(BitsIn,varargin)
if nargin==1 % Works as an interleaver
map=randperm(length(BitsIn));
varargout{1}=BitsIn(map);
varargout{2}=map;
else % Works as a deinterleaver
map=varargin{1};
BitsOut=BitsIn; 
BitsOut(map)=BitsIn;  
varargout{1}=BitsOut;
end
end