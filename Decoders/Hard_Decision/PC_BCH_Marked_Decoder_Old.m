%% PC decoder with miscorrection detection based on zero-syndrome and marked bits 
function hat_coded_bits=PC_BCH_Marked_Decoder(coded_bits,BCHnu,BCHt,BCHe,It,llr,R,varargin)
% Decoding function for a product code.
% It uses the function eBCHdecoder (from Yi Lei) which assumes the BCH code
% was extended by one bit. BCH codes with no extension or with more than
% one extension bits (C. Haeger uses 2 for example) has not been
% implemented yet.
% M is a mask (logical array) of marked bits (based on their reliability) for each received
% block

if nargin>7
ItThresh=varargin{1};                      % Iterations threshold for zero-syndrome decoding rejection
end

BCHn=(2^BCHnu)-1;                % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;     % Information bits
dmin=2*BCHt+BCHe+1;                % for t=2, e=1

    % Sort and mark bits based on reliability threshold
    [~,Ir]=sort(abs(llr).');
    Ir=Ir.';
    [~,Ic]=sort(abs(llr));
    MarkedBits=(abs(llr)>R);   % Marks bits based on LLRs


hat_coded_bits=coded_bits;
%CheckConflict=zeros(1,BCHn+BCHe);
ZeroSyndromeIdx=zeros(1,BCHn+BCHe); % Initialise as none of the codewords having 0 syndrome
status=zeros(1,BCHn+BCHe);

for ll=1:It
    
    % Decode Rows
    for rr=1:BCHn+BCHe
        MiscDetFlag=0;  
        Fail=0;
        
        [~, tmp_decoded, status(rr)]=eBCHdecoder(hat_coded_bits(rr,:),BCHnu,BCHk,BCHt);
        
       if ~(status(rr)==0 || status(rr)==1)     % Decoding failure
        Fail=1;  % Flag for failure event
        hat_coded_bits(rr,Ir(rr,1))=~hat_coded_bits(rr,Ir(rr,1));    % Flips least reliable bit   
        [~, tmp_decoded, status(rr)]=eBCHdecoder(hat_coded_bits(rr,:),BCHnu,BCHk,BCHt); % New decoding attempt 
        

        if ~(status(rr)==0 || status(rr)==1)  % Additional failure after bit flip
          hat_coded_bits(rr,Ir(rr,1))=~hat_coded_bits(rr,Ir(rr,1));    % Flips back least reliable bit     
          continue
        end
       end 
        
        % Miscorrection detection for first ItThresh+1 iterations 
       if ll<ItThresh+1
        while MiscDetFlag<2    
        Conflict=(ZeroSyndromeIdx & ~(tmp_decoded==hat_coded_bits(rr,:))) | (MarkedBits(rr,:) & ~(tmp_decoded==hat_coded_bits(rr,:)));
        CheckConflict=~isempty(find(Conflict,1));
        
        if CheckConflict  % detects miscorrection
           MiscDetFlag= MiscDetFlag+1;
            % detects miscorrection
            % Bit flipping
           if ~Fail 
               
               if MiscDetFlag<2  % Flip bits after miscorrection detected
                we=sum(tmp_decoded~=hat_coded_bits(rr,:));  % Error weight
                Nflip=dmin-we-BCHt;
                hat_coded_bits(rr,(Ir(rr,1:Nflip)))=~hat_coded_bits(rr,(Ir(rr,1:Nflip))); % Bit flipping
                [~, tmp_decoded, status(rr)]=eBCHdecoder(hat_coded_bits(rr,:),BCHnu,BCHk,BCHt); % Second decoding attempt     
                
               else 
                hat_coded_bits(rr,(Ir(rr,1:Nflip)))=~hat_coded_bits(rr,(Ir(rr,1:Nflip))); % Flips bits back if miscorrection is still detected 
                status(rr)=2;
                
               end
           else
           
           hat_coded_bits(rr,Ir(rr,1))=~hat_coded_bits(rr,Ir(rr,1));    % Flips back least reliable bit     
           break %  exit miscorrection detection loop 
           end
           
           
        else  % no miscorrection
        hat_coded_bits(rr,:)=tmp_decoded;    
        %status(rr)=0;  
        break % exit miscorrection detection loop 
        end
        
        end % end of while loop
        
        else % succesful decoding
        hat_coded_bits(rr,:)=tmp_decoded;
        %status(rr)=0; 
        end  % end of miscorrection detection iteration loop
       
    end
    
     ZeroSyndromeIdx=(status==0)|(status==1);
     %ZeroSyndromeIdx=(status==0);

    % Decode Columns
    for cc=1:BCHn+BCHe 
        MiscDetFlag=0;  
        Fail=0;
        [~, tmp_decoded, status(cc)]=eBCHdecoder(hat_coded_bits(:,cc).',BCHnu,BCHk,BCHt);
        
       
      if ~(status(cc)==0 || status(cc)==1)     % Decoding failure
        
        hat_coded_bits((Ic(1,cc)),cc)=~hat_coded_bits((Ic(1,cc)),cc);    % Flips least reliable bit   
        [~, tmp_decoded, status(cc)]=eBCHdecoder(hat_coded_bits(:,cc).',BCHnu,BCHk,BCHt);
        Fail=1;  % Flag for failure event
        
        if ~(status(cc)==0 || status(cc)==1)
          hat_coded_bits((Ic(1,cc)),cc)=~hat_coded_bits((Ic(1,cc)),cc);    % Flips back least reliable bit     
          continue 
        end
        
      end
        
      if ll<ItThresh                 % Miscorrection detection for first ItThresh-1 iterations        
               
        while MiscDetFlag<2  
        Conflict=(ZeroSyndromeIdx & ~(tmp_decoded==hat_coded_bits(:,cc))) | (MarkedBits(:,cc) & ~(tmp_decoded==hat_coded_bits(:,cc)));
        CheckConflict=~isempty(find(Conflict,1));
        
        if CheckConflict  % Decode if no conflict
         MiscDetFlag= MiscDetFlag+1;
            % detects miscorrection
            % Bit flipping
            if ~Fail 
               
             if MiscDetFlag<2  % Flip bits after miscorrection detected
       
                we=sum(tmp_decoded~=hat_coded_bits(:,cc));  % Error weight
                Nflip=dmin-we-BCHt;
                hat_coded_bits((Ic(1:Nflip,cc)),cc)=~hat_coded_bits((Ic(1:Nflip,cc)),cc); % Bit flipping
             [~, tmp_decoded, status(cc)]=eBCHdecoder(hat_coded_bits(:,cc).',BCHnu,BCHk,BCHt); % Second decoding attempt            
             else 
                hat_coded_bits((Ic(1:Nflip,cc)),cc)=~hat_coded_bits((Ic(1:Nflip,cc)),cc); % Flips bits back if miscorrection is still detected 
                status(cc)=2;
            end
            else  % miscorrection with previous decoding failure
           hat_coded_bits((Ic(1,cc)),cc)=~hat_coded_bits((Ic(1,cc)),cc);    % Flips back least reliable bit     
           break 
           end
        else  % succesful decoding
        hat_coded_bits(:,cc)=tmp_decoded.';    
        %status(cc)=0;  
        break
        end
        
        end % end of while loop
        
      else 
       hat_coded_bits(:,cc)=tmp_decoded.';
        
      end
        
    end
    
     ZeroSyndromeIdx=(status==0)|(status==1);
     %ZeroSyndromeIdx=(status==0);

  %  hat_coded_bits(Mark)=coded_bits(Mark);    % keep reliable bits as they were before the decoding at the end of each iteration
end