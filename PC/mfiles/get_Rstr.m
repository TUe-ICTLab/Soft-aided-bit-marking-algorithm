function R_str   = get_Rstr(R)

switch R 
    case 1/4, R_str='1_4';
    case 1/3, R_str='1_3';
    case 2/5, R_str='2_5';  %0.4
    case .42, R_str='0_42'; %0.42
    case .45, R_str='0_45'; %0.45
    case .47, R_str='0_47'; %0.47
	case 1/2, R_str='1_2';  %0.5
    case .53, R_str='0_53'; %0.53
    case .57, R_str='0_57'; %0.57
	case 3/5, R_str='3_5';  %0.6
    case .62, R_str='0_62'; %0.62
    case .64, R_str='0_64'; %0.64
    case 2/3, R_str='2_3';  %0.66
    case .69, R_str='0_69'; %0.69
    case .71, R_str='0_71'; %0.71
    case .72, R_str='0_72'; %0.72
	case 3/4, R_str='3_4';
    case .77, R_str='0_77'; %0.77
    case 4/5, R_str='4_5';
    case .81, R_str='0_81'; %0.81
    case .82, R_str='0_82'; %0.82
    case 5/6, R_str='5_6';  %0.83
    case .84, R_str='0_84';
    case .86, R_str='0_86';
    case .87, R_str='0_87';
    case 8/9, R_str='8_9';  %0.88
    case 9/10,R_str='9_10'; %0.90
    case .91, R_str='0_91';
    case .92, R_str='0_92';
    case .94, R_str='0_94';
end
