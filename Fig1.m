

sx = 0.7;
mx = 1.4;
x = exp(randn(100,1)*sx+mx);
y = poissrnd(20*x);


u = randn(100,1);

su1 = sx/2;
su2 = sx;

x1 = x.*exp(u*su1-1/2*su1^2);
x2 = x.*exp(u*su2-1/2*su2^2);
xi=linspace(0,max(x2),1000);


figure(2);clf
subplot(131)
plot(x,y/10,'o','markersize',3)
hold all
xlim([0 max(x2)])
refline(2,0)
ylabel({'recruitment','(seedlings, saplings)'})
title('no error')

subplot(132)
plot(x1,y/10,'o','markersize',3)
hold all
phat = glmfit(log(x1),y,'poisson');
plot(xi,(exp(phat(1))*xi.^phat(2))/10,'r-')
xlim([0 max(x2)])
refline(2,0)
xlabel('conspecific density (adults or seeds)')
title('small error')
legend('observed data','fitted line','true line');legend('boxoff')

subplot(133)
plot(x2,y/10,'o','markersize',3)
hold all
phat = glmfit(log(x2),y,'poisson');
plot(xi,(exp(phat(1))*xi.^phat(2))/10,'r-')
xlim([0 max(x2)])
refline(2,0)
title('large error')

sublabel
