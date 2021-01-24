function hat_coded_bits=PC_BCH_Decoder(coded_bits,BCHnu,BCHt,BCHe,L)
% Decoding function for a product code.
% It uses the function eBCHdecoder (from Yi Lei) which assumes the BCH code
% was extended by one bit. BCH codes with no extension or with more than
% one extension bits (C. Haeger uses 2 for example) has not been
% implemented yet.

BCHn=(2^BCHnu)-1;                % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;     % Information bits

if (BCHe==0 || BCHe>2),return;end % Not implemented (
hat_coded_bits=coded_bits;
for ll=1:L
    % Decode Rows
    for rr=1:BCHn+BCHe
        [tmp_err_poly, tmp_decoded, status]=eBCHdecoder(hat_coded_bits(rr,:),BCHnu,BCHk,BCHt);
        hat_coded_bits(rr,:)=tmp_decoded;
    end
    % Decode Columns
    for cc=1:BCHn+BCHe
        [tmp_err_poly, tmp_decoded, status]=eBCHdecoder(hat_coded_bits(:,cc).',BCHnu,BCHk,BCHt);
        hat_coded_bits(:,cc)=tmp_decoded.';
    end
end