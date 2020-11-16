%=========================HW3 problem C (c)==============================
clear;
k = 7;
N = 2^k;
delta_t = 4/N^2;
timepoint = N^2/128;
U = zeros(N+1,N+1,timepoint+1);

tic


[D,x] = cheb(N);
I = eye(N-1);
D2 = D^2;
D2 =D2(2:N,2:N);
L = kron(I,D2)+kron(D2,I);
%clear I D2
L_square = L^2;
f = @(x)...
        exp(-400*(x)^2);
lap_f = @(x)...
    640000*x^2*exp(-400*x^2) - 800*exp(-400*x^2);
Ut = zeros(N+1,N+1);
lap_Ut = zeros(N+1,N+1);
%======================fine grid=========================================
for i = 1:N+1
    for j = 1:N+1
        Ut(i,j) = f(x(i))*f(x(j));
        lap_Ut(i,j) = lap_f(x(i))*f(x(j))+f(x(i))*lap_f(x(j)); 
    end
end
U(:,:,2) = delta_t*(Ut+delta_t^2/6*lap_Ut);

for n = 3:timepoint+1
    U(2:N,2:N,n) = 2*U(2:N,2:N,n-1)-U(2:N,2:N,n-2)+delta_t^2*...
        reshape(L*reshape(U(2:N,2:N,n-1),[],1),[],N-1) + ...
        1/12*delta_t^4*...
        reshape(L_square*reshape(U(2:N,2:N,n-1),[],1),[],N-1); %correction term
end
    
clear L L_square
error_vec = zeros(3,1);
%==========================comparing with fine grid======================
for k = 4:6
    N = 2^k;
    delta_t = 4/N^2;
    timepoint = N^2/128;
    Unew = zeros(N+1,N+1,timepoint);
    [D,x] = cheb(N);
    I = eye(N-1);
    D2 = D^2;
    D2 =D2(2:N,2:N);
    L = kron(I,D2)+kron(D2,I);

    L_square = L^2;
    Ut = zeros(N+1,N+1);
    lap_Ut = zeros(N+1,N+1);
    
    for i = 1:N+1
        for j = 1:N+1
            Ut(i,j) = f(x(i))*f(x(j));
            lap_Ut(i,j) = lap_f(x(i))*f(x(j))+f(x(i))*lap_f(x(j)); 
        end
    end
    Unew(:,:,2) = delta_t*(Ut+delta_t^2/6*lap_Ut);

    for n = 3:timepoint+1
        Unew(2:N,2:N,n) = 2*Unew(2:N,2:N,n-1)-Unew(2:N,2:N,n-2)+delta_t^2*...
            reshape(L*reshape(Unew(2:N,2:N,n-1),[],1),[],N-1) + ...
            1/12*delta_t^4*...
            reshape(L_square*reshape(Unew(2:N,2:N,n-1),[],1),[],N-1);%correction term
    end
    error_vec(k-3) = norm(reshape(Unew-U(1:2^(7-k):129,...
        1:2^(7-k):129,1:4^(7-k):129),[],1),inf);
end
    
plot((4:6),-log2(error_vec),'*')
line((4:6),-log2(error_vec),'LineStyle','-')
line([4,6],[16,24]);


toc






%===========================function Chebb=================================
function [D,x] = cheb(N)
if N == 0
    D= 0; 
    x=1; 
    return, 
end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D - diag(sum(D,2));
% off-diagonal entries
end

