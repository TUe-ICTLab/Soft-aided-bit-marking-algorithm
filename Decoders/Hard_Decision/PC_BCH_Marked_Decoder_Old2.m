function hat_coded_bits=PC_BCH_Marked_Decoder2(coded_bits,BCHnu,BCHt,BCHe,L,varargin)
% Decoding function for a product code.
% It uses the function eBCHdecoder (from Yi Lei) which assumes the BCH code
% was extended by one bit. BCH codes with no extension or with more than
% one extension bits (C. Haeger uses 2 for example) has not been
% implemented yet.
% M is a mask (logical array) of marked bits (based on their reliability) for each received
% block

llr=varargin{1};   % Block of LLRs passed to the decoder
R=varargin{2};     % LLR threshold for reliability
%flip=varargin{3};   % Flag to enable bit-flipping

BCHn=(2^BCHnu)-1;                % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;     % Information bits
ItThresh=2;                      % Iterations threshold for zero-syndrome decoding rejection

if ~isempty(varargin{1})
    % Sort and mark bits based on reliability threshold
    [~,Ir]=sort(abs(llr.'));
    [~,Ic]=sort(abs(llr));
    MarkedBits=(abs(llr)>R);   % Marks bits based on LLRs
end

hat_coded_bits=coded_bits;
%CheckConflict=zeros(1,BCHn+BCHe);
ZeroSyndromeIdx=zeros(1,BCHn+BCHe); % Initialise as none of the codewords having 0 syndrome
status=zeros(1,BCHn+BCHe);
%flag0=1;
%flag1=0;

for ll=1:L    
    % Decode Rows
    for rr=1:BCHn+BCHe
        [tmp_err_poly, tmp_decoded, status(rr)]=eBCHdecoder(hat_coded_bits(rr,:),BCHnu,BCHk,BCHt);
        
        % Decoding failure
        if ~(isequal(status(rr),1) || isequal(status(rr),0)) % Flip most unreliable bit in case of 
           hat_coded_bits(rr,Ir(1,rr))=~hat_coded_bits(rr,Ir(1,rr));
        
        % Miscorrection detection for first ItThresh+1 iterations 
        elseif ll<=ItThresh
        %elseif ll<-1
        %if ll==1
        Conflict=(ZeroSyndromeIdx & ~(tmp_decoded==hat_coded_bits(rr,:))) | (MarkedBits(rr,:) & ~(tmp_decoded==hat_coded_bits(rr,:)));
        %else 
        %Conflict=(ZeroSyndromeIdx & ~(tmp_decoded==hat_coded_bits(rr,:)));
        %end
        
        CheckConflict=~isempty(find(Conflict,1));
        
        if CheckConflict  % detects miscorrection
       
            status(rr)=2;
           % flips the first d0-w(e)-t most unreliable bits 
%             if flip && ll==1
%             disp('Boh');
%             we=sum(tmp_err_poly);  % Error weight
%             flick=BCHt+2-we;
%             hat_coded_bits(rr,Ir(1:flick,rr))=~hat_coded_bits(rr,Ir(1:flick,rr));
%             end
%          

        else   % Successful decoding with no miscorrection detected
        hat_coded_bits(rr,:)=tmp_decoded;    
        end
             
        else % Successful decoding
        hat_coded_bits(rr,:)=tmp_decoded;
        end
     end   
      %flag0=1;   
      ZeroSyndromeIdx=((status==0)|(status==1));
     
     %ZeroSyndromeIdx=(status==0);

    % Decode Columns
    for cc=1:BCHn+BCHe 
        [tmp_err_poly, tmp_decoded, status(cc)]=eBCHdecoder(hat_coded_bits(:,cc).',BCHnu,BCHk,BCHt);
        
        if ~(isequal(status(cc),1) || isequal(status(cc),0)) 
           hat_coded_bits(Ic(1,cc),cc)=~hat_coded_bits(Ic(1,cc),cc); % Flip most unreliable bit in case of decoding failure
        
        
        % Miscorrection detection enabled for first ItThresh iterations 
        elseif ll<ItThresh 
        %if ll==1    
        Conflict=(ZeroSyndromeIdx & ~(tmp_decoded==hat_coded_bits(:,cc))) | (MarkedBits(:,cc) & ~(tmp_decoded==hat_coded_bits(:,cc)));
        %else
        %Conflict=(ZeroSyndromeIdx & ~(tmp_decoded==hat_coded_bits(:,cc)));
        %end

        CheckConflict=~isempty(find(Conflict,1));
        
        %% Miscorrection detection
        if CheckConflict   
           
           status(cc)=2;          % Miscorrection detected
%            if flip && ll==1         % Flip unreliable bits
%            we=sum(tmp_err_poly);  % Error weight
%            flick=BCHt+2-we;
%            hat_coded_bits(Ic(1:flick,cc),cc)=~hat_coded_bits(Ic(1:flick,cc),cc);
%            end
%            
%            % Set flags for bit flipping 
%            if flag0   
%            flag1=flag1+1;
%            flag0=0;
%            end
%            
        else % Successful decoding with no miscorrection detected
        hat_coded_bits(:,cc)=tmp_decoded.';       
           
        end
        
        
        else % Successful decoding
        hat_coded_bits(:,cc)=tmp_decoded.';    
        
        end
        
        %else 
        %hat_coded_bits(:,cc)=tmp_decoded.';
           
     end
     %flag0=1;  
     ZeroSyndromeIdx=((status==0)|(status==1));
     %ZeroSyndromeIdx=(status==0);

  %  hat_coded_bits(Mark)=coded_bits(Mark);    % keep reliable bits as they were before the decoding at the end of each iteration
end

end

