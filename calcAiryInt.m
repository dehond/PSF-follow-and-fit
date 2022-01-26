%calcAiryInt.m
%Calculates the volume underneath the Airy disk given NA and wavelength.
function integralValue = calcAiryInt(NA, lambda)
    k = 2*pi/lambda;
    dr = 0.01;
    r = dr:dr:100;
    airy = (2*besselj(1, k*NA*r)./(k*NA*r)).^2;

    integralValue = 2*pi*trapz(r, r.*airy);
end