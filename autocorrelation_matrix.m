function [ Ra ] = autocorrelation_matrix( Ra_temp )
%AUTOCORRELATION_MATRIX Summary of this function goes here
%   Detailed explanation goes here

%         autocorrelation matrix Ra
Ra_acf   = xcorr(Ra_temp);
Rxx_temp = Ra_acf(length(Ra_temp):end);
Ra       = toeplitz(Rxx_temp,[Rxx_temp(1) conj(Rxx_temp(2:end))]);


end

