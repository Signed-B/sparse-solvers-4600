clc; clear; close all;

n=10; %number points per direction
x1=0;
x2=1;
y1=0;
y2=1;

function [A,rhs,D]=poisson_setup(x1,x2,y1,y2,n,ddx,ddy)
    %A is matrix, rhs is randhand side, D is domain of finite difference
    %points
    %x1,x2,y1,y2 define interval of interest
    %n is number of point in each direction
    %ddx,ddy are double partials in each direction

    %assumes function is zero at end conditions
    A=zeros(n^2,n^2);
    A=A++4*eye(n^2,n^2);
    rhs=zeros(n^2,1);
    for i=1:n^2
        
    end


end




function [z]=f(x,y)
    z=sin(2*pi*x)+sin(2*pi*y);
end

function [z]=ddx(x,y) %d''f/dx''
    z=-4*pi^2*sin(2*pi*x);
end


function [z]=ddy(x,y) %d''f/dy''
    z=-4*pi^2*sin(2*pi*y);
end
