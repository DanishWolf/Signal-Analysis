function [Type] = sigtest(s, Samp)
% if length(s) ~= Samp
%     error("Please chek samples, singla may not be complete")
% end

if mod(length(s),Samp) ~= 0
    error("Please chek samples, singal may not be complete")
end
% threshold of arbitrary samll number
thresh = 1e-12;

%Centres the function to correct calculation
s = s(2:end);


if s- flip(s) < thresh
    Type = "Even";
    
elseif s + flip(s) < thresh
    Type = "Odd";
    
else
    Type = "Neither";
end

end

