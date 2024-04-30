clc; clear; close all;

n=11; %number points per direction
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

    Dx=linspace(x1,x2,n);
    Dy=linspace(y1,y2,n);
    D=zeros(n^2,2);
    D(1,1)=x1;
    for i=1:n^2-1
        D(i+1,1)=Dx(mod(i,n)+1);
        D(i+1,2)=Dy(floor(i/n)+1);
    end
    

    
    rhs=zeros(n^2,1);
    rhs=D
    for i=1:n^2
        if i<n
            rhs(i)=0;
        end
        
        if i>n^2-n
            rhs(i)=0;
        end
    end


        A=zeros(n^2,n^2);


end




function [z]=f(x,y) %solution
    z=sin(2*pi*x)+sin(2*pi*y);
end

function [z]=ddx(x,y) %d''f/dx''
    z=-4*pi^2*sin(2*pi*x);
end


function [z]=ddy(x,y) %d''f/dy''
    z=-4*pi^2*sin(2*pi*y);
end
