function [ output ] = myMapping( c,mod )
%MYMAPPING Summary of this function goes here
%   Detailed explanation goes here

switch mod
    case 2 % QPSK
        % Mapping 
        % 00 -> 1+1i
        % 01 -> +1-1i
        % 10 -> -1+1i
        % 11 -> -1-1i
        mapping_vector = [1+1i +1-1i -1+1i -1-1i]./sqrt(2);
        j = 1;
        for i = 1:2:length(c)
            if c(i) == 0 && c(i+1) == 0
                output(j) = mapping_vector(1);
            elseif c(i) == 0 && c(i+1) == 1
                output(j) = mapping_vector(2);
            elseif c(i) == 1 && c(i+1) == 0
                output(j) = mapping_vector(3);
            else
                output(j) = mapping_vector(4);
            end
            j = j + 1;
        end
    
end

end

