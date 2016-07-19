num_of_iter = 1e5;

cond_D0 = zeros(1,num_of_iter);

a = 4; b = 5;
J_all = a + (b-a).*rand(1,num_of_iter);

a = 1e-3; b = 0.25;
J_per_N_all = a + (b-a).*rand(1,num_of_iter);

N_all = round(J_all./J_per_N_all);

index = find((N_all.*J_per_N_all)<4);
N_all(index) = N_all(index) + 1;

index = find((N_all.*J_per_N_all)>5);
N_all(index) = N_all(index) - 1;

a = 0; b = 20000;
A_all = a + (b-a).*rand(1,num_of_iter);

a = 0; b = 20000;
B_all = a + (b-a).*rand(1,num_of_iter);



for i=1:num_of_iter;

J = J_all(i);

J_per_N = J_per_N_all(i);

N = N_all(i);

A = A_all(i);

B = B_all(i);


D0 =  [cos(2*pi*J_per_N*(0:(N-1)))'  sin(2*pi*J_per_N*(0:(N-1)))'  ones(N,1)];

D0_4 = -A*(0:(N-1)).*sin(2*pi*J_per_N.*(0:(N-1))) + B*(0:(N-1)).* cos(2*pi*J_per_N.*(0:(N-1)));

scale = sqrt(A^2 + B^2)*N/1.852;

D0 = [D0 D0_4'/scale];

cond_D0(i) = cond(D0'*D0/N);

end

subplot(1,2,1)
 hist(cond_D0,12:0.01:17)
xlim([12 17])
ylim([0 600])
xlabel ({'Condition number of ({\bfD}{\it_i}^T{\bfD}{\it_i}){\it_s_c}';'(a)'})
ylabel ('Number of occurances')

cond_D0 = zeros(1,num_of_iter);

a = 100; b = 101;
J_all = a + (b-a).*rand(1,num_of_iter);

a = 1e-3; b = 0.25;
J_per_N_all = a + (b-a).*rand(1,num_of_iter);

N_all = round(J_all./J_per_N_all);

index = find((N_all.*J_per_N_all)<4);
N_all(index) = N_all(index) + 1;

index = find((N_all.*J_per_N_all)>5);
N_all(index) = N_all(index) - 1;

a = 0; b = 20000;
A_all = a + (b-a).*rand(1,num_of_iter);

a = 0; b = 20000;
B_all = a + (b-a).*rand(1,num_of_iter);



for i=1:num_of_iter;

J = J_all(i);

J_per_N = J_per_N_all(i);

N = N_all(i);

A = A_all(i);

B = B_all(i);


D0 =  [cos(2*pi*J_per_N*(0:(N-1)))'  sin(2*pi*J_per_N*(0:(N-1)))'  ones(N,1)];

D0_4 = -A*(0:(N-1)).*sin(2*pi*J_per_N.*(0:(N-1))) + B*(0:(N-1)).* cos(2*pi*J_per_N.*(0:(N-1)));

scale = sqrt(A^2 + B^2)*N/1.852;

D0 = [D0 D0_4'/scale];

cond_D0(i) = cond(D0'*D0/N);

end

subplot(1,2,2)
 hist(cond_D0,12:0.01:17)
xlim([12 17])
ylim([0 13000])
xlabel ({'Condition number of ({\bfD}{\it_i}^T{\bfD}{\it_i}){\it_s_c}';'(b)'})
ylabel ('Number of occurances')



