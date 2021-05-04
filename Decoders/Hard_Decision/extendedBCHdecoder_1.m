%% BDD-extended
function output1=extendedBCHdecoder_1(enc,dec,s,t,input1)

output1=zeros(1,length(input1));

[cv,error_idicate]=step(dec,[zeros(s,1);input1(1,1:end-1)']);

if error_idicate==-1
    
 output1=input1;
 
else
     
 de1=sum(input1)+error_idicate;
 de=mod(de1,2);
 
 if (de+error_idicate)>t
  
 output1=input1;    
     
 else
 
 cvp=step(enc,cv);    
 output1(1:end-1)=cvp(s+1:end)';   
 output1(end)=mod(input1(end)+de,2);    
     
 end
    
 
end

end