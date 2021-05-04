%% iBDD-SR decoder (Alireza.s implementation)
%%% Example: decodeout1=y_decode_func(vec,receivechan,256,239,8,2,8,0,coeffvec,enc,dec)
%%% vec: Product code array
%%% receivechan: Channel reliability

function decodeout1=iBDD_SR(vec,receivechan,block_comp,inf_comp,itermax1,itermax2,v,s,coeffvec,enc,dec)

old_vec=vec;
new_vec=vec;

channelout1=(receivechan);
channelout2=(receivechan).';

t_code=(block_comp-inf_comp-1)/v;


ind_vec=1;

for iter=1:itermax1
   
coeff=coeffvec(ind_vec);    
ind_vec=ind_vec+1;

for i=1:block_comp
    
[indx,out1]=extendedBCHdecoder(enc,dec,s,t_code,vec(i,:));

if (indx)               % bounded distance decoding for Reed solomon codes
    
aaa=(sign(channelout1(i,:))+1)/2;
vec(i,:)=aaa; 
    
else
    
bbb=2*(out1)-1;
aaa=(sign(coeff*bbb+channelout1(i,:))+1)/2;
vec(i,:)=aaa;


end

end

vec2=vec';

% coeff=coeffvec(ind_vec); 
% ind_vec=ind_vec+1;

for i=1:block_comp
    
[indx,out2]=extendedBCHdecoder(enc,dec,s,t_code,vec2(i,:));

if (indx)               % bounded distance decoding for Reed solomon codes
    
aaa=(sign(channelout2(i,:))+1)/2;
vec2(i,:)=aaa; 


else
   
bbb=2*(out2)-1;
aaa=(sign(coeff*bbb+channelout2(i,:))+1)/2;    
vec2(i,:)=aaa;

end

end

vec=vec2';
new_vec=vec;

if length(find(new_vec~=old_vec))==0
    
    break;
    
else
    
old_vec=vec;

end

end


%% iBDD-phase
for iter=1:itermax2


for i=1:block_comp
    
out1=extendedBCHdecoder_1(enc,dec,s,t_code,vec(i,:));    
vec(i,:)=out1;   

end

vec2=vec';

for i=1:block_comp
    
out2=extendedBCHdecoder_1(enc,dec,s,t_code,vec2(i,:));
vec2(i,:)=out2;

end

vec=vec2';

new_vec=vec;

if length(find(new_vec~=old_vec))==0
    
    break;
    
else
    
old_vec=vec;

end

end

decodeout1=vec;

end
%% 
