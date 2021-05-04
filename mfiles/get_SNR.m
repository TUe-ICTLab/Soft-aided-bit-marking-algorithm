function EsN0dB=get_SNR(m,R)
% This function returns precalculated values of SNR for different
% modulation formats and rates
% This is adjusted "manually".
%
% Alex Alvarado
% June 2017

if m==1
    modulation='PM-QPSK',     
else
   modulation=['PM-',num2str(2^(2*m)),'QAM'];
end
    
switch R
    case 1/4
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-7.5:.1:-5.4];
            case 'PM-16QAM',    EsN0dB=[-3.5:.1:-1.2];
            case 'PM-64QAM',    EsN0dB=[-0.5:.1:2.1];
            case 'PM-256QAM',   EsN0dB=[2.5:.1:4.8];
            case 'PM-1024QAM',  EsN0dB=[4:.1:7];
        end
    case 1/3
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-5.4:.1:-4.3];
            case 'PM-16QAM',    EsN0dB=[-1.5:.1:0.4];
            case 'PM-64QAM',    EsN0dB=[[2:.1:3.9],[3.95]];
            case 'PM-256QAM',   EsN0dB=[5:.1:7.2];
        end
    case 2/5
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[-4.3:.1:-3.5],[-3.45]];
            case 'PM-16QAM',    EsN0dB=[0:.1:1.5];
            case 'PM-64QAM',    EsN0dB=[3.5:.1:5.6];
            case 'PM-256QAM',   EsN0dB=[7.4:.1:9.2];
        end
	case 0.42
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[-4:.1:-3.15],[-3.2]];
        end
	case 0.45
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-3.8:.1:-2.7];
        end
	case 0.47
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-3.5:.1:-2.4];
        end
    case 1/2
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-2.9:.1:-2.1];
            case 'PM-16QAM',    EsN0dB=[1.5:.1:3.2];
            case 'PM-64QAM',    EsN0dB=[5:.1:7.6];
            case 'PM-256QAM',   EsN0dB=[10:.1:11.8];
        end
	case 0.53
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-2.5:.1:-1.6];
        end
    case 0.57
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-2:.1:-1.2];
        end
	case 3/5
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[-1.6:.1:-0.9],[-0.85]];
            case 'PM-16QAM',  	EsN0dB=[3.1:.1:4.8]; 
            case 'PM-64QAM',    EsN0dB=[7.5:.1:9.5];
            case 'PM-256QAM',   EsN0dB=[[12.6:.1:14.0]];
        end
	case 0.62
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[-1.5:.1:-0.6],[-0.55]];
        end
	case 0.64
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-1.4:.1:-0.3];
        end
     case 2/3
         switch modulation
            case 'PM-QPSK',     EsN0dB=[[-1.2:.1:-0.1],[-0.15,-0.02,0]];
            case 'PM-16QAM',    EsN0dB=[[4.5:.1:5.7],[5.72,5.74]]; 
            case 'PM-64QAM',    EsN0dB=[9:.1:10.6];
         end
	case 0.69
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-.8:.1:0.3];
            case 'PM-16QAM',    EsN0dB=[4.6:.1:6.1];
        end
    case 0.71
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-.6:.1:0.6];
            case 'PM-16QAM',    EsN0dB=[4.7:.1:6.5];
            case 'PM-64QAM',    EsN0dB=[10.5:.1:11.6];
        end
    case 0.72
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-.5:.1:0.8];
            case 'PM-16QAM',    EsN0dB=[4.8:.1:6.7];
        end
     case 3/4
         switch modulation
            case 'PM-QPSK',     EsN0dB=[[0:.1:0.9],[0.95,0.97]];
            case 'PM-16QAM',    EsN0dB=[5:.1:7];
            case 'PM-64QAM',    EsN0dB=[9.5:.1:12.2];
            case 'PM-256QAM',   EsN0dB=[15.5:.1:17.5];
         end
    case 0.77
        switch modulation
            case 'PM-QPSK',     EsN0dB=[-.4:.1:1.3];
            case 'PM-16QAM',    EsN0dB=[5.5:.1:7.4];
        end
     case 4/5
         switch modulation
             case 'PM-QPSK',	EsN0dB=[[1.0:.1:1.5],[1.55,1.57,1.59]];
             case 'PM-16QAM',   EsN0dB=[[6.0:.1:7.7],[7.72,7.74]];
%             case 'PM-64QAM',    EsN0dB=[10.5:.1:13.1];
         end
    case 0.81
        switch modulation
            case 'PM-QPSK',     EsN0dB=[0.5:.1:1.8];
            case 'PM-16QAM',    EsN0dB=[5.5:.1:8];
            case 'PM-64QAM',    EsN0dB=[11:.1:13.4];
        end
    case 0.82
        switch modulation
            case 'PM-QPSK',     EsN0dB=[1.2:.1:2.0];
            case 'PM-16QAM',    EsN0dB=[6.8:.1:8.2];
            case 'PM-64QAM',    EsN0dB=[12.8:.1:13.6];
        end
    case 5/6
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[1.5:.1:2.0],[2.05,2.09]];
            case 'PM-16QAM',    EsN0dB=[[7.0:.1:8.3],[8.32,8.34,8.36]];
            case 'PM-64QAM',    EsN0dB=[12.5:.1:13.8];
            case 'PM-256QAM',   EsN0dB=[18.2:.1:19.4];
        end
    case 0.84
        switch modulation
            case 'PM-QPSK',     EsN0dB=[1.8:.1:2.3];
            case 'PM-16QAM',    EsN0dB=[7.6:.1:8.6];
        end
    case 0.86
        switch modulation
            case 'PM-QPSK',     EsN0dB=[1:.1:2.6];
            case 'PM-16QAM',    EsN0dB=[7:.1:8.8];
            case 'PM-64QAM',    EsN0dB=[12:.1:14.5];
        end
    case 0.87
        switch modulation
            case 'PM-QPSK',     EsN0dB=[2.0:.1:2.8];
            case 'PM-16QAM',    EsN0dB=[7.8:.1:9.2];
            case 'PM-64QAM',    EsN0dB=[13.7:.1:14.7];
        end
     case 8/9
         switch modulation
             case 'PM-QPSK',   	EsN0dB=[[2:.1:3],[3.05,3.1,3.13,3.15]];
             case 'PM-16QAM', 	EsN0dB=[[8:.1:9.4],[9.46,9.5,9.52]];
%             case 'PM-64QAM',    EsN0dB=[13.5:.1:15.1];
         end
    case 9/10
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[2:.1:3.3],[3.34]];
            case 'PM-16QAM',    EsN0dB=[[8:.1:9.7],[9.72,9.74,9.76]];
            case 'PM-64QAM',    EsN0dB=[[13.5:.1:15.4],[15.35,15.42,15.44,15.45,15.46]];
            case 'PM-256QAM',   EsN0dB=[19.7:.1:21.1];
        end
    case 0.91
        switch modulation
            case 'PM-QPSK',     EsN0dB=[2.6:.1:3.5];
            case 'PM-16QAM',    EsN0dB=[8.8:.1:10];
            case 'PM-64QAM',    EsN0dB=[15:.1:15.7];
        end
    case 0.92
        switch modulation
            case 'PM-QPSK',     EsN0dB=[2.7:.1:3.9];
            case 'PM-16QAM',    EsN0dB=[9.2:.1:10.4];
        end
    case 0.94
        switch modulation
            case 'PM-QPSK',     EsN0dB=[[2.8:.1:4.3],[4.35]];
            case 'PM-16QAM',    EsN0dB=[9.8:.1:10.8];
        end
end