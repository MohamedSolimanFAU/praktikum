function filename = generate_channel(channel_type,N_r,N_t,Anz_channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HISTORY:
% 26.08.08: Channels for LTE (EPA, EVA, ETU), VA30 added
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename_prefix = 'LTE_channel_' ;
save_type    = 'cell' ; 'vector' ;

switch channel_type
    % E-UTRA propagation conditions
    case 'EPA' %Extended Pedestrian A
        pdp_t_ns    = [0 30 70 90 110 190 410];
        pdp         = 10.^([0 -1 -2 -3 -8 -17.2 -20.8]/10);        
    case 'EVA' %Extended Vehicular A
        pdp_t_ns    = [0 30 150 310 370 710 1090 1730 2510];
        pdp         = 10.^([0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.7 -12.0 -16.9]/10);
    case 'ETU' %Extended Typical Urban
        pdp_t_ns    = [0 50 120 200 230 500 1600 2300 5000];
        pdp         = 10.^([-1.0 -1.0 -1.0 0 0 0 -3 -5 -7]/10);
    % Equalizer test
    case 'EQ10'
        pdp_t_ns    = (0:1:10-1)*100 ;       
        pdp         = ones(size(pdp_t_ns)) ;             % Varianzen 
    case 'EQ100'
        pdp_t_ns    = (0:1:100-1)*100 ;       
        pdp         = ones(size(pdp_t_ns)) ;             % Varianzen
        
    case 'EQ144'
        pdp_t_ns    = (0:1:144-1)*100 ;       
        pdp         = ones(size(pdp_t_ns)) ;             % Varianze  
    
    % Lilly's original channel DA (Diplom Arbeit)
    case 'DA'
        pdp         = [ 1/2 1/3 1/6 ] ;             % Varianzen
        pdp_t_ns    = [   0 100 200 ] ;
   
    % UTRA propagation conditions
    % ITU Pedestrian A Speed 3km/h (PA3)
    case 'A'
        pdp         = 10.^( [  0 -9.7 -19.2 -22.8 ] /10) ;
        pdp_t_ns    = [  0  100  200 400 ] ;
    % ITU Pedestrian B Speed 3km/h (PB3)
    case 'B'
        pdp         = 10.^( [  0 -0.9 -4.9 -8.0 -7.8 ] /10) ;
        pdp_t_ns    = [  0  200  800 1200 2300 ] ;
    
    case 'VA30' % Vehicular A 30km/h
        pdp_t_ns    = [ 0  310 710 1090 1730 2510 ] ;
        pdp         = 10.^( [0 -1 -9 -10 -15 -20]/10) ;
   
    case 'F'
        pdp         = 1 ;
        pdp_t_ns    = 0 ;
        
    case 'exp'
        zeit_norm = [0:1:144];
        zeit_0=-144/log(0.1);
        pdp = exp(-1*zeit_norm/zeit_0);
        pdp_t_ns    = 100*zeit_norm;
        
    case 'C_lang'
        pdp         = 10.^( [  0 -1.0 -9.0 -10.0 -15.0 -20.0 ] /10) ;
        pdp_t_ns    = [  0 300 700 1100 1700 2500 ] ;
        
    case 'D'
        pdp         = 10.^( [  -2.5 0 -12.8 -10.0 ] /10) ;
        pdp_t_ns    = [  0 300 8900 12900 ] ;
        
end ;



filename      = [ filename_prefix channel_type '_Anz' num2str(Anz_channel) '_' save_type '_NR' num2str(N_r) '_NT' num2str(N_t) '.mat' ] ;

pdp           = pdp/ sum(pdp) ;                     % Normierung
%Anz_pdp       = length(pdp) ;
q_h           = max(ceil(pdp_t_ns/100)) ;

pos_w         = round(pdp_t_ns/100)+1 ;
w_h           = zeros(1,q_h+1) ;
w_h(pos_w)    = sqrt(pdp) ;                         % Gewichte 

channel = zeros(N_r,N_r,q_h+1);
h_vek = zeros(N_r,N_t,Anz_channel,q_h+1);
%g_vek = zeros(1,1,Anz_channel,q_h+1);
for i = 1:N_r
    for l = 1:N_t
        g_vek                           = ( randn(Anz_channel,q_h+1) + j*randn(Anz_channel,q_h+1) ) / sqrt(2) ;
        h_vek(i,l,:,1:q_h+1)            = g_vek .* repmat(w_h,Anz_channel,1) ;
        temp                            = h_vek(i,l,:,:);
        avg_power                       = sum(abs(h_vek(:)).^2)/Anz_channel ;
        channel                         = h_vek(i,l,:,:) / sqrt(avg_power) ;           % Normierung 
    end
end
switch save_type
    case 'vector'    
        save(filename,'h_vek','h_vek','q_h','q_h','Anz_channel','Anz_channel') ;
    case 'cell'
        h_cell        = cell(1,Anz_channel) ;
        for i = 1:Anz_channel
            h_cell{i} = reshape(h_vek(:,:,i,:),N_r,N_t,q_h+1) ;            
        end ; 
       save(filename,'h_cell','h_cell','q_h','q_h','Anz_channel','Anz_channel','pdp','pdp_t_ns') ;        
end ;


%% Aufruf
%
% %INIT
% load(filename) ;
% 
% %LOOP akt_channel_no = 1:Anz_channel
% h = h_cell{akt_channel_no}
