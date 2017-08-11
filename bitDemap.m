%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soft-Demapper for QPSK und 16QAM (AWGN-channel)                             %
% Based on the diploma thesis of Uyen Ly Dang                                 %
% File creation date: 07.04.2008                                              %
% Author            : Michael Ruder                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function[output, LLR]= softDemapping(u, mod, noiseVar, llr_fig)
% 23.07.2012: File wurde umbenannt, sonst nix (DANG)
function[output, LLR]= bitDemap(u, mod, noiseVar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:        mod:      modulation scheme: QPSK or 16 QAM                   %
%               u:        transmitted (and equalized)symbols (row vector)     %
%               noiseVar: estimated noise variance                            %
%               (llr_fig:  set llr_fig=1 to plot the LLR)                     %
% OUTPUT:       output:   demapped data - interleaved coded bits output       %
%               LLR:      corresponding  LLR (row vector)                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mod = '16QAM';
%noiseVar = 0.01;
llr_fig=0;
%u = u.';
switch mod
    case 2 % == QPSK
        %u_Im = sign(randn(13,1))
        %u_Re = sign(randn(13,1))
        %u = u_Re+j*u_Im;

        % Initialization
        LLR = zeros(1,2);
        c0 = zeros(2,2);
        c1 = zeros(2,2);

        % c0(mu,:) Vektor of symbols where c[mu] == 0
        c0 =  [1+1i 1-1i;
                    1+1i -1+1i]./sqrt(2);
        c1 = [-1+1i -1-1i;
              1-1i -1-1i]./sqrt(2);
        
        %c0(2,:) = [1+j -1+j]./sqrt(2);
        %c1(2,:) = [1-j -1-j]./sqrt(2);

        % Plot LLR of QPSK
        if llr_fig==1
            LLR = zeros(1,1,2);
            %LLLR = zeros(1,1,2);
            im = -1.3:0.01:1.3;
            re = -1.3:0.01:1.3;
            for ind_r = 1:length(re)
                for ind_i = 1:length(im)
                    condPdf_c0 = zeros(1,2);
                    condPdf_c1 = zeros(1,2);
                    for temp = 1:2
                        condPdf_c0(1) = condPdf_c0(1) + exp(-(abs(complex(re(ind_r),im(ind_i))-c0(1,temp))).^2./noiseVar);
                        condPdf_c1(1) = condPdf_c1(1) + exp(-(abs(complex(re(ind_r),im(ind_i))-c1(1,temp))).^2./noiseVar);
                        condPdf_c0(2) = condPdf_c0(2) + exp(-(abs(complex(re(ind_r),im(ind_i))-c0(2,temp))).^2./noiseVar);
                        condPdf_c1(2) = condPdf_c1(2) + exp(-(abs(complex(re(ind_r),im(ind_i))-c1(2,temp))).^2./noiseVar);
                    end
                    LLR(ind_r,ind_i,1:2)=log(condPdf_c1./condPdf_c0);
                end
            end
            figure; mesh(re,im,LLR(:,:,1));title('LLR[1]');xlabel('Re'); ylabel('Im'), zlabel('LLR[1]');
            figure; mesh(re,im,LLR(:,:,2));title('LLR[2]');xlabel('Re'); ylabel('Im'), zlabel('LLR[2]');
            clear re; clear im; clear ind_r; clear ind_i;
            LLR = zeros(1,2);
            condPdf_c0 = zeros(1,2);
            condPdf_c1 = zeros(1,2);
        end


        % Compute LLR        
        U = zeros(2*length(u),2);
        
        % Old version by Lilly
        %for k = 1:length(u)
        %    U(2*k-1:2*k,1:2) = ones(2,2).*u(k);
        %end
        % New version by Michael Ruder
        temp = [u u; u u];
        U    = reshape(temp, length(u)*2, 2);
        
        condPdf_c0 = realmin + sum( exp(-abs(U-repmat(c0,length(u),1)).^2./noiseVar), 2).';
        condPdf_c1 = realmin + sum( exp(-abs(U-repmat(c1,length(u),1)).^2./noiseVar), 2).';
        
        LLR = reshape(log(condPdf_c1 ./ condPdf_c0), 2, length(u)).';
        
        c_temp = sign(LLR(1:length(u),:).');
%       idx = find(c_temp < 0);
        c_temp(c_temp<0) = 0;
        LLR = LLR.'; LLR = LLR(:).'; %% <-- ??! 
        output(1:2*length(u)) = c_temp(:);
        
    case '3' % ==8PSK
    case '4' % == 16QAM
        %u_Im = [sign(randn(13,1)); 3*sign(randn(13,1))];
        %idx = randperm(length(u_Im));
        %u_Im = u_Im(idx);
        %u_Re = [sign(randn(13,1)); 3*sign(randn(13,1))];
        %idx = randperm(length(u_Im));
        %u_Re = u_Re(idx);
        % u = u_Re+j*u_Im;

        %Initialization
        LLR = zeros(1,4);
        % c0(mu,:) Vektor of symbols where c[mu] == 0
        c0 = [1+3*j 1+j 1-j 1-3*j  3+3*j  3+j   3-j   3-3*j ;
             -3+j  -1+j 1+j 3+j   -3+3*j -1+3*j 1+3*j 3+3*j ;
              1+3*j 1+j 1-j 1-3*j -1+3*j -1+j  -1-j  -1-3*j ; 
             -3+j  -1+j 1+j 3+j   -3-j   -1-j   1-j   3-j ] ./sqrt(10);
        
        c1 = [-1+3*j -1+j  -1-j  -1-3*j -3+3*j -3+j  -3-j  -3-3*j ;
              -3-j   -1-j   1-j   3-j   -3-3*j -1-3*j 1-3*j 3-3*j ;
               3+3*j  3+j   3-j   3-3*j -3+3*j -3+j  -3-j  -3-3*j ;
              -3+3*j -1+3*j 1+3*j 3+3*j -3-3*j -1-3*j 1-3*j 3-3*j ] ./sqrt(10);

        %Plot LLR
        if llr_fig==1
            im = -3.3:0.01:3.3;
            re = -3.3:0.01:3.3;
            for ind_r = 1:length(re)
                for ind_i = 1:length(im)
                    condPdf_c0 = zeros(1,4);
                    condPdf_c1 = zeros(1,4);
                    for temp = 1:8
                        %--------------------------------------------------------%
                        % Das noch mal vielleicht kompakter
                        % ------------------------------------------------
                        condPdf_c0(1) = condPdf_c0(1) + exp(-(abs(complex(re(ind_r),im(ind_i))-c0(1,temp))).^2./noiseVar);
                        condPdf_c1(1) = condPdf_c1(1) + exp(-(abs(complex(re(ind_r),im(ind_i))-c1(1,temp))).^2./noiseVar);
                        condPdf_c0(2) = condPdf_c0(2) + exp(-(abs(complex(re(ind_r),im(ind_i))-c0(2,temp))).^2./noiseVar);
                        condPdf_c1(2) = condPdf_c1(2) + exp(-(abs(complex(re(ind_r),im(ind_i))-c1(2,temp))).^2./noiseVar);
                        condPdf_c0(3) = condPdf_c0(3) + exp(-(abs(complex(re(ind_r),im(ind_i))-c0(3,temp))).^2./noiseVar);
                        condPdf_c1(3) = condPdf_c1(3) + exp(-(abs(complex(re(ind_r),im(ind_i))-c1(3,temp))).^2./noiseVar);
                        condPdf_c0(4) = condPdf_c0(4) + exp(-(abs(complex(re(ind_r),im(ind_i))-c0(4,temp))).^2./noiseVar);
                        condPdf_c1(4) = condPdf_c1(4) + exp(-(abs(complex(re(ind_r),im(ind_i))-c1(4,temp))).^2./noiseVar);

                    end
                    ind0 = find(condPdf_c0 ==0); condPdf_c0(ind0)=0.000000000001;
                    ind1 = find(condPdf_c1 ==0); condPdf_c1(ind1)=0.000000000001;
                    clear ind0;
                    clear ind1;
                    LLR(ind_r,ind_i,1:4)=log(condPdf_c1./condPdf_c0);
                end
            end
            figure;mesh(re,im,LLR(:,:,1));title('LLR[1]');xlabel('Re'); ylabel('Im'), zlabel('LLR[1]');
            figure;mesh(re,im,LLR(:,:,2));title('LLR[2]');xlabel('Re'); ylabel('Im'), zlabel('LLR[2]');
            figure;mesh(re,im,LLR(:,:,3));title('LLR[3]');xlabel('Re'); ylabel('Im'), zlabel('LLR[3]');
            figure;mesh(re,im,LLR(:,:,4));title('LLR[4]');xlabel('Re'); ylabel('Im'), zlabel('LLR[4]');
            clear re; clear im; clear ind_r; clear ind_i;
            LLR = zeros(1:4);
            condPdf_c0 = zeros(1,4);
            condPdf_c1 = zeros(1,4);
        end

        % Compute LLR
        U = zeros(4*length(u),8);

        % Old version by Lilly
        %for k = 1:length(u)
        %    U(4*k-3:4*k,1:8) = ones(4,8).*u(k);
        %end
        % New Version by Michael Ruder
        temp = [u u u u u u u u; u u u u u u u u; u u u u u u u u; u u u u u u u u ];
        U    = reshape(temp, length(u)*4, 8);
        
        condPdf_c0 = realmin + sum( exp(-abs(U-repmat(c0,length(u),1)).^2./noiseVar), 2).';
        condPdf_c1 = realmin + sum( exp(-abs(U-repmat(c1,length(u),1)).^2./noiseVar), 2).';
        
        LLR = reshape(log(condPdf_c1 ./ condPdf_c0), 4, length(u)).';

        c_temp = sign(LLR(1:length(u),:).');
%       idx = find(c_temp < 0);
        c_temp(c_temp<0) = 0;
        LLR = LLR.'; LLR = LLR(:).'; %% <-- ??! 
        output(1:4*length(u)) = c_temp(:);        
end