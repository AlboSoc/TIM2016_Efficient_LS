sigma1 = 0.5;
sigma2 = 1;

z_all = 2:0.001:5;
cond_num = zeros(1,length(z_all));

for i = 1:length(z_all)

z = z_all(i);
    
sigma3 = (z/6 + 0.5 + sqrt(z^2/36 + z/12 + 0.25))/2;

sigma4 = (z/6 + 0.5 - sqrt(z^2/36 + z/12 + 0.25))/2;

sigma_max = max([sigma1 sigma2 sigma3 sigma4]);

sigma_min = min([sigma1 sigma2 sigma3 sigma4]);

cond_num(i) = sigma_max/sigma_min;

end

plot(z_all,cond_num,'LineWidth',2)
ylim([13.9 19])
ylabel ('Condition number of ({\bfD}{\it_i}^T{\bfD}{\it_i}){\it_s_c}')
xlabel ('{\itz} = {\itR}^2{\itN}^2/\gamma^2')