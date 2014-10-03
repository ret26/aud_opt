function YHWR = softHWR(Y,rho)


YHWR = log(1+exp(Y*rho))/rho; % half wave rectify
indLg = Y*rho>100;
YHWR(indLg) = Y(indLg);
