clc; clear; close all;

n=6; %number points per direction

x1=0; %bounrdries
x2=1;
y1=0;
y2=1; 

%assumes constant space in both direction (with break otherwise)
h=(x2-x1)/(n-1); %spacing

%Domain (100 different ways)
[grid_x,grid_y]=meshgrid(x1:h:x2,y1:h:y2);
Dx=linspace(x1,x2,n);
Dy=linspace(y1,y2,n);
D=zeros(n^2,2);
D(1,1)=x1; 
for i=1:n^2-1
    D(i+1,1)=Dx(mod(i,n)+1);
    D(i+1,2)=Dy(floor(i/n)+1);
end


A=zeros(n^2,n^2); %creates empty variables
rhs=zeros(n^2,1);

for i=1:n^2
    on_boundry=0; %boundry boolean
    if i<=n || i>n^2-n || mod(i,n)==0 || mod(i,n)==1  %checks if u_i is on boundry
        on_boundry=1;
        rhs(i)=g(D(i,1),D(i,2)); %if yes rhs is boundry function
        A(i,i)=1; %corresponding row in A is just part of identity
    else
        rhs(i)=f(D(i,1),D(i,2)); %if not rhs is f
        for j=1:n^2
            A(i,i)=-4/h^2; %approximates second derivvative
            A(i,i+1)=1/h^2;
            A(i,i-1)=1/h^2;
            A(i,i+(n))=1/h^2;
            A(i,i-(n))=1/h^2;
        end
    end
end

U=A\rhs; %solves system

grid_u=zeros(n,n); %converst solution back to matrix form
for i=1:n
    for j=1:n
        grid_u(i,j)=U(n*(i-1)+j);
    end
end
hold on;
surface(grid_x,grid_y,grid_u) %plots solution

% calcs and plot actual
[real_x,real_y]=meshgrid(x1:(x2-x1)/200:x2,y1:(x2-x1)/200:y2);
surface(real_x,real_y,actual(real_x,real_y))


function [z]=f(x,y) %rhs forcing function
    z=-8*pi^2*sin(2*pi*y).*sin(2*pi*x);
end

function [z]=g(x,y) %broundry condition function
    z=0; % zero for dirichlet boundry conditions
end

function [z]=actual(x,y) %actual solution to ivp (for plots and error)
    z=sin(2*pi*x).*sin(2*pi*y);
end