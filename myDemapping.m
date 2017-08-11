function [ output ] = myDemapping( c,mod )
%MYDEMAPPING Summary of this function goes here
%   Detailed explanation goes here

switch mod
    case 2 % QPSK
        % Mapping 
        % 00 -> 1+1i
        % 01 -> +1-1i
        % 10 -> -1+1i
        % 11 -> -1-1i
        j = 1;
        for i = 1:length(c)
            if real(c(i)) >= 0 && imag(c(i)) >= 0
                output(j:j+1) = [0 0];
            elseif real(c(i)) >= 0 && imag(c(i)) < 0
                output(j:j+1) = [0 1];
            elseif real(c(i)) < 0 && imag(c(i)) >= 0
                output(j:j+1) = [1 0];
            else
                output(j:j+1) = [1 1];
            end
            j = j + 2;
        end
end

