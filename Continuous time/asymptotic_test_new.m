%% asymptotic expansion test

% set parameters
pi1 = 0.6;
pi2 = 1-pi1;
u = 1e3; % k*tau, assumed to be very large

% set functions
tol = 1e-8;
wVec = linspace(tol,1-tol,1e5); % w in (0,1)

% P(w,s_n=2|s_{n-1}=1,theta)
Our_func = @(w) pi1*sqrt(u/(4*pi))*1./(pi1*pi2*w.*(1-w)).^(1/4).*exp(-u*(sqrt(pi2*w)-sqrt(pi1*(1-w))).^2);

% P(w,s_n=2|s_{n-1}=1,theta)
Our_func2 = @(w) pi1*sqrt(u/(4*pi))*sqrt(w./(1-w)).*1./(pi1*pi2*w.*(1-w)).^(1/4).*exp(-u*(sqrt(pi2*w)-sqrt(pi1*(1-w))).^2);

% exact formula for P(w,s_n=2|s_{n-1}=1,theta)
I0_1_func = @(w) pi2*u*exp(-u*(pi2*w+pi1*(1-w))).*besseli(0,2*u*sqrt(pi1*pi2*w.*(1-w)));

% exact formula for P(w,s_n=2|s_{n-1}=1,theta)
I0_2_func = @(w) pi1*u*exp(-u*(pi2*w+pi1*(1-w))).*besseli(0,2*u*sqrt(pi1*pi2*w.*(1-w)));

% exact formula for P(w,s_n=1|s_{n-1}=1,theta)
I1_1_func = @(w) u*exp(-u*(pi2*w+pi1*(1-w))).*sqrt(pi1*pi2*w./(1-w)).*besseli(1,2*u*sqrt(pi1*pi2*w.*(1-w)));

% exact formula for P(w,s_n=2|s_{n-1}=2,theta)
I1_2_func = @(w) u*exp(-u*(pi2*w+pi1*(1-w))).*sqrt(pi1*pi2*(1-w)./w).*besseli(1,2*u*sqrt(pi1*pi2*w.*(1-w)));

% P(w|theta) formula (2.12) in Berezhkovskii with pre-factor
Berez_func = @(w) sqrt(u/(4*pi))*1./(pi1*pi2*w.*(1-w)).^(1/4).*exp(-u*(sqrt(pi2*w)-sqrt(pi1*(1-w))).^2);

% gaussian distribution with pre-factor
gauss_func = @(w,mu,sigma) 1./sqrt(2*pi*sigma^2).*exp(-(w-mu).^2/(2*sigma^2));

% test similarity to a gaussian
mu = pi1;
sigma = sqrt(2*pi1*pi2/u);

yGauss = gauss_func(wVec,mu,sigma);
% trapz(wVec,yGauss) % test normalisation
% yOur = Our_func(wVec);
yOur2 = Our_func2(wVec);
yI0_1 = I0_1_func(wVec);
yI0_2 = I0_2_func(wVec);
yI1_1 = I1_1_func(wVec);
yI1_2 = I1_2_func(wVec);
yBerez = Berez_func(wVec);

% plot(wVec,yOur2,'r-','LineWidth',1.5)
plot(wVec,yI0+yI1_2,'r-','LineWidth',1.5)
hold on
plot(wVec,yGauss,'b--','LineWidth',1.5)
hold on
% plot(wVec,yBerez,'k-','LineWidth',1.5)





