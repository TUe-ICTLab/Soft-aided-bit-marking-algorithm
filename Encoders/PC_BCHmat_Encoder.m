function coded_bits=PC_BCHmat_Encoder(info_bits,enc,BCHnu,BCHt,BCHe)
% This function implements a produc code encoder. The component code is a
% BCH code with opne bit extension. No extension or more than one bit
% extension (C. Haeger uses some of those codes) has not been implemented
% yet (but it should be trivial).
%
% info_bits should be a (2^nu)-1-nu*t by (2^nu)-1-nu*t vector of bits
% coded_bits will be a 1 by (2^nu)-1+e vector of bits

BCHn=(2^BCHnu)-1;                           % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;                % Information bits
coded_bits=zeros(BCHn+BCHe,BCHn+BCHe);      % Matrix of coded bits

if BCHe>2,return;end % Not implemented
%% Row Encoding
ext_bits=[];
for rr=1:BCHk
    %tmp_coded_bits=bchenc_mex(BCHnu,BCHn,BCHt,info_bits(rr,:)); % Encode
    tmp_coded_bits=step(enc,info_bits(rr,:).').';
    if BCHe==1 % sum mod2 all the 2^nu-1 bits
        ext_bits=mod(sum(tmp_coded_bits),2);
    elseif BCHe==2 % sum mod2 all even and odd bits
        ext_bits=[mod(sum(tmp_coded_bits(2:2:end-BCHe)),2),mod(sum(tmp_coded_bits(1:2:end-BCHe)),2)];
    end
    coded_bits(rr,:)=[tmp_coded_bits,ext_bits]; % Append extension bits
end
%% Column Encoding, including parities on parities
new_info_bits=coded_bits(1:BCHk,:); % copy the first k rows
for cc=1:BCHn+BCHe
     %tmp_coded_bits=bchenc_mex(BCHnu,BCHn,BCHt,new_info_bits(:,cc).'); % Encode
     tmp_coded_bits=step(enc,new_info_bits(:,cc)).';

    if BCHe==1 % sum mod2 all the 2^nu-1 bits
        ext_bits=mod(sum(tmp_coded_bits),2);
    elseif BCHe==2 % sum mod2 all even and odd bits
        ext_bits=[mod(sum(tmp_coded_bits(2:2:end-BCHe)),2),mod(sum(tmp_coded_bits(1:2:end-BCHe)),2)];
    end
    coded_bits(:,cc)=[tmp_coded_bits,ext_bits].'; % Append extension bits (here we rewrite the info bits. not needed but who cares...)
end
return

