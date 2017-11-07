%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for generating and saving channel realisations of LTE MIMO channels  %
% Based on the modified script by Partick Nickel                              %
% File creation date: 03.04.2008                                              %
% Author            : Michael Ruder                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filename = generate_channel(channel_type,N_r,N_t,number_of_channel_realisations)
filename_prefix = 'LTE_channel_' ;
%N_r          = 2;
%N_t          = 4;
%channel_type = 'A' ;'B' ; 
%number_of_channel_realisations  = 1000 ;
save_type    = 'cell' ; 'vector' ;

switch channel_type
        % inverse linearly decaying pdp
    case 'linInv'
        l_h = 5;
        l = 1:l_h;
        pdp         = 1./l ;             % Varianzen
        pdp_t_ns    = (0:l_h-1).*100 ;   
    % This case is for Equalizer Testing 10 Tabs
    case 'EQ10'
        pdp_t_ns    = (0:1:10-1)*100 ;       
        pdp         = ones(size(pdp_t_ns)) ;             % Varianzen 
    % This case is for Equalizer Testing 100 Tabs
    case 'EQ100'
        pdp_t_ns    = (0:1:100-1)*100 ;       
        pdp         = ones(size(pdp_t_ns)) ;             % Varianzen 
    % The following case 'A' is not pedestrian A --> Lilly's original channel DA=  Diplom Arbeit
    case 'Lilly-DA'
        pdp         = [ 1/2 1/3 1/6 ] ;             % Varianzen
        pdp_t_ns    = [   0 100 200 ] ;
    % This case is ITU Pedestrian A Speed 3km/h (PA3)
    case 'ITU-PA'
        pdp         = 10.^( [  0 -9.7 -19.2 -22.8 ] /10) ;
        pdp_t_ns    = [  0  100  200 400 ] ;
    % This case is ITU Pedestrian B Speed 3km/h (PB3)
    case 'ITU-PB'
        pdp         = 10.^( [  0 -0.9 -4.9 -8.0 -7.8 ] /10) ;
        pdp_t_ns    = [  0  200  800 1200 2300 ] ;
    % This case is AWGN Channel --> Flat Fading    
    case 'Flat-Fading'
        pdp         = 1 ;
        pdp_t_ns    = 0 ;
    case 'ETU'
        pdp         = 10.^( [  -1 -1.3 0 0 0 -3 -5 -7 ] /10) ;
        pdp_t_ns    = [  0  100  200 300 500 1600 2300 5000 ] ;
    otherwise
    error('ERROR: The channel type %c is not supported!\n', channel_type);
end ;



filename      = [ filename_prefix channel_type '_Anz' num2str(number_of_channel_realisations) '_' save_type '_NR' num2str(N_r) '_NT' num2str(N_t) '.mat' ] ;
pdp           = pdp/ sum(pdp) ;                     % Normierung
%Anz_pdp       = length(pdp) ;
q_h           = max(ceil(pdp_t_ns/100)) ;

pos_w         = round(pdp_t_ns/100)+1 ;
w_h           = zeros(1,q_h+1) ;
w_h(pos_w)    = sqrt(pdp) ;                         % Gewichte 

channel = zeros(N_r,N_r,q_h+1);
h_vek = zeros(N_r,N_t,number_of_channel_realisations,q_h+1);
%g_vek = zeros(1,1,number_of_channel_realisations,q_h+1);
for i = 1:N_r
    for l = 1:N_t
        g_vek                           = ( randn(number_of_channel_realisations,q_h+1) + 1i*randn(number_of_channel_realisations,q_h+1) ) / sqrt(2) ;
        h_vek(i,l,:,1:q_h+1)            = g_vek .* repmat(w_h,number_of_channel_realisations,1) ;
        temp                            = h_vek(i,l,:,:);
        avg_power                       = sum(abs(h_vek(:)).^2)/number_of_channel_realisations ;
        channel                         = h_vek(i,l,:,:) / sqrt(avg_power) ;           % Normierung 
    end
end
switch save_type
    case 'vector'    
        save(filename,'h_vek','h_vek','q_h','q_h','number_of_channel_realisations','number_of_channel_realisations') ;
    case 'cell'
        h_cell        = cell(1,number_of_channel_realisations) ;
        for i = 1:number_of_channel_realisations
            h_cell{i} = reshape(h_vek(:,:,i,:),N_r,N_t,q_h+1) ;            
        end ; 
        save(filename,'h_cell','q_h','number_of_channel_realisations') ;        
end ;


%% How to use this:
%
% %INIT
% load(filename) ;
% 
% %LOOP akt_channel_no = 1:number_of_channel_realisations
% h = h_cell{akt_channel_no}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HISTORY:                                                                    %
% 03.04.2008 : File created.                                                  %
% 21.04.2008 : Add number_of_channel_realistions to save file name.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
