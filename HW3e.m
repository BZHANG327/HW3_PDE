%=============================HW3 problem C (e)==========================
clear;
tic
k = 6;
N = 2^k;
delta_t = 12/N^2;
timepoint = N^2/8;
U = zeros(N+1,N+1,timepoint+1);
Utrue = zeros(N+1,N+1,timepoint+1);
[D,x] = cheb(N);
I = eye(N-1);
D2 = D^2;
D2 =D2(2:N,2:N);
L = kron(I,D2)+kron(D2,I);
%clear I D2
L_square = L^2;
error_vec = zeros(5);
%==========================FD scheme=======================================
for B = [2,4,8,16,32]
    U_x = @(x,B)...
        sin(B*pi*x);
    U_t = @(x,B)...
        sin(sqrt(2)*pi*B*x);
    deltax = 1/N;
    for i = 1:N+1
        for j = 1:N+1
            for k = 1:timepoint+1
                Utrue(i,j,k) = U_x((i-1)*deltax,B)*U_x((j-1)*deltax,B)*U_t((k-1)*delta_t,B)/...
                    (sqrt(2)*B*pi);
            end
        end
    end
    for i = 1:N+1
        for j = 1:N+1
            U(i,j,2) = delta_t*U_x((i-1)*deltax,B)*U_x((j-1)*deltax,B);
        end
    end
    for n = 3:timepoint+1
        for i = 2:N
            for j = 2:N
                U(i,j,n) = (delta_t^2/deltax^2)*...
                    (U(i+1,j,n-1)+U(i,j+1,n-1)+U(i-1,j,n-1)+U(i,j-1,n-1)-4*U(i,j,n-1))...
                    -U(i,j,n-2)+2*U(i,j,n-1);
            end
        end
    end
    error_vec(round(log2(B))) = norm(reshape(U-Utrue,[],1),inf);
end
hold on 
plot([2,4,8,16,32],log10(error_vec),'*')
line([2,4,8,16,32],log10(error_vec),'LineStyle','-')
%==============================
for B = [2,4,8,16,32]
    U_x = @(x,B)...
        sin(B*pi*x);
    U_t = @(x,B)...
        sin(sqrt(2)*pi*B*x);
    for i = 1:N+1
        for j = 1:N+1
            for k = 1:timepoint+1
                Utrue(i,j,k) = U_x(x(i),B)*U_x(x(j),B)*U_t((k-1)*delta_t,B)/...
                    (sqrt(2)*B*pi);
            end
        end
    end
    Ut = zeros(N+1,N+1);
    lap_Ut = zeros(N+1,N+1);
    for i = 1:N+1
        for j = 1:N+1
            Ut(i,j) = U_x(x(i),B)*U_x(x(j),B);
            lap_Ut(i,j) = -2*B^2*pi^2*Ut(i,j); 
        end
    end
    U(:,:,2) = delta_t*(Ut+delta_t^2/6*lap_Ut);
    for n = 3:timepoint+1
        U(2:N,2:N,n) = 2*U(2:N,2:N,n-1)-U(2:N,2:N,n-2)+delta_t^2*...
            reshape(L*reshape(U(2:N,2:N,n-1),[],1),[],N-1) + ...
            1/12*delta_t^4*...
            reshape(L_square*reshape(U(2:N,2:N,n-1),[],1),[],N-1); %correction term
    end
    error_vec(round(log2(B))) = norm(reshape(U-Utrue,[],1),inf);
end
plot([2,4,8,16,32],log10(error_vec),'o')
line([2,4,8,16,32],log10(error_vec),'LineStyle','-')
line([0,35],[-3,-3]);
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