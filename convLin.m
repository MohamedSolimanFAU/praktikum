%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear convolution of MIMO channels.                                        %
% File creation date: 03.12.2008                                              %
% Author            : Michael Ruder                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% output = convLin_mr(input, channel)
% input    : [N_t x blocklength] data sequence
% channel  : [N_r x N_t x q_h+1] MIMO channel in the time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output   : [N_r x blocklength + q_h] convolved MIMO signal
%            Note: "real" data starts with the q_h+1 - th symbol!
% 

function output = convLin(input,... % Datensequenz [N_t x length(input)]
                            channel)  % Übertragungskanal [N_r x N_t x q_h+1]
                      
% Linear Convolution

N_r       = size(channel,1);
N_t       = size(channel,2);
q_h       = size(channel,3)-1;
l_input   = size(input,2);
output    = zeros(N_r,l_input+q_h);

for N_R = 1:N_r
  for N_T = 1:N_t
    output(N_R,:) = output(N_R,:) + conv(input(N_T,:), shiftdim(channel(N_R, N_T, :),1));
  end
end