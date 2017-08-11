% Gray-Mapping and Modulation for 4QAM and 16QAM

function[output,mappingVector] = bitMap(c,mod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:     
% c: interleaved coded bits (row_vector)
% mod: modulation scheme QPSK or 16QAM, 64QAM
% 
% OUTPUT:    
% output:    sequence of complex symbols (row vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mod = '16QAM';
%idx = 0;

%c = [zeros(1,7) ones(1,8)];
%idx = randperm(length(c));
%c = c(idx);

switch mod
    case 1 % BPSK
        mappingVector = 0;
        p = 0;
        output = 0;
        c_t = zeros(1,ceil(length(c)));
        c_t(1,1:length(c)) = c;
        
        
        mappingVector = [1 -1];
        display('ab hier keine BPSK-Unterstützung')
    case 2 % bits per symbol (QPSK)
        % Initialization
        mappingVector = 0;
        p = 0;
        output = 0;
        
        % "Zero-Padding"
        c_t = zeros(1,ceil(length(c)/2)*2);
        c_t(1,1:length(c)) = c; 
        
        % Mapping 
        % 00 -> 1+1i
        % 01 -> +1-1i
        % 10 -> -1+1i
        % 11 -> -1-1i
        mappingVector = [ 1+1i +1-1i -1+1i -1-1i]./sqrt(2);
        output = ones(size(c_t, 1), length(c_t)/2);
        for k = 1:2           
            output = output + 2^(2-k) * c_t(:, k:2:length(c_t));            
        end
        output = mappingVector(output);
    case 3 % 8 PSK        
        c_t = zeros(1,ceil(length(c)));
        c_t(1,1:length(c)) = c;
        output = ones(size(c_t, 1), length(c_t)/3);
        % bin to dec mapping
        for k = 1:3           
            output = output + 2^(3-k) * c_t(:, k:3:length(c_t));            
        end
        
        mappingVector = [exp(1i*2*pi*0/8) exp(1i*2*pi*1/8) exp(1i*2*pi*3/8) exp(1i*2*pi*2/8) exp(1i*2*pi*7/8) exp(1i*2*pi*6/8) exp(1i*2*pi*4/8) exp(1i*2*pi*5/8)];        
        output = mappingVector(output);
    case 4 % 16QAM
        % Initialization
        mappingVector = 0;
        p = 0;
        output = 0;
        
        % "Zero-Padding"
        c_t = zeros(1,ceil(length(c)/4)*4);
        c_t(1,1:length(c)) = c;

        %Mapping
        mappingVector = [1+1i 1+3*1i 3+1i 3+3*1i 1-1i 1-3*1i 3-1i 3-3*1i -1+1i -1+3*1i -3+1i -3+3*1i -1-1i -1-3*1i -3-1i -3-3*1i]./ sqrt(10);                
        
        output = ones(size(c_t, 1), length(c_t)/4);
        %% convert from bits to symbol numbers
        
        for k=1:4           
            output = output + 2^(4-k) * c_t(:, k:4:length(c_t));
        end
        output = mappingVector(output) ;
    case 6
       % Initialization
        mappingVector = 0;
        p = 0;
        output = 0;
        
        % "Zero-Padding"
        c_t = zeros(1,ceil(length(c)/6)*6);
        c_t(1,1:length(c)) = c; 
        
        
        %Mapping
        mappingVector = [3+3*1i  3+1*1i  1+3*1i  1+1*1i  3+5*1i  3+7*1i  1+5*1i  1+7*1i  5+3*1i  5+1*1i  7+3*1i  7+1*1i  5+5*1i  5+7*1i  7+5*1i  7+7*1i,... % bit i 
                         3-3*1i  3-1*1i  1-3*1i  1-1*1i  3-5*1i  3-7*1i  1-5*1i  1-7*1i  5-3*1i  5-1*1i  7-3*1i  7-1*1i  5-5*1i  5-7*1i  7-5*1i  7-7*1i,...
                        -3+3*1i -3+1*1i -1+3*1i -1+1*1i -3+5*1i -3+7*1i -1+5*1i -1+7*1i -5+3*1i -5+1*1i -7+3*1i -7+1*1i -5+5*1i -5+7*1i -7+5*1i -7+7*1i,... % bit i 
                        -3-3*1i -3-1*1i -1-3*1i -1-1*1i -3-5*1i -3-7*1i -1-5*1i -1-7*1i -5-3*1i -5-1*1i -7-3*1i -7-1*1i -5-5*1i -5-7*1i -7-5*1i -7-7*1i]./ sqrt(42);
        
        output = ones(size(c_t, 1), length(c_t)/4);
        %% convert from bits to symbol numbers
        
        for k=1:6           
            output = output + 2^(6-k) * c_t(:, k:6:length(c_t));
        end
        output = mappingVector(output) ;
end


%EOF

%{
SymbTab_32QAM       = [-0.6708 - 1.1180i  -0.2236 - 1.1180i  -0.6708 + 1.1180i  -0.2236 + 1.1180i  -1.1180 - 0.6708i    -1.1180 - 0.2236i  -1.1180 + 0.6708i  -1.1180 + 0.2236i  -0.2236 - 0.6708i  -0.2236 - 0.2236i     -0.2236 + 0.6708i  -0.2236 + 0.2236i  -0.6708 - 0.6708i  -0.6708 - 0.2236i  -0.6708 + 0.6708i     -0.6708 + 0.2236i   0.6708 - 1.1180i   0.2236 - 1.1180i   0.6708 + 1.1180i   0.2236 + 1.1180i     1.1180 - 0.6708i   1.1180 - 0.2236i   1.1180 + 0.6708i   1.1180 + 0.2236i   0.2236 - 0.6708i     0.2236 - 0.2236i   0.2236 + 0.6708i   0.2236 + 0.2236i   0.6708 - 0.6708i   0.6708 - 0.2236i    0.6708 + 0.6708i   0.6708 + 0.2236i];
SymbTab_8PSK        = exp(2*pi*j*[3 4 2 1 6 5 7 0]/8) ;
SymbTab_BPSK        = exp(2*pi*j*[0 1]/2) ;
SymbTab_QPSK        = [1+j +1-j -1+j -1-j];
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HISTORY
% 29.04.2008 Ausgabe des Mapping Vektors
%