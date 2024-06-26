clc; clear; close all; %Housekeeping
%%
n=5; %number points per direction

x1=0; %bounrdries must be sqaure
x2=1;
y1=0;
y2=1; 

%function handles for forcing, end condition, and actual solution
F=@(x,y)f(x,y);
G=@(x,y)g(x,y);
Actual=@(x,y)actual(x,y);

[A,rhs,grid_x,grid_y]=create_system_P2D(n,x1,x2,y1,y2,G,F); %creates sparse banded matrix to solve problem

U=A\rhs; %solves system (default matlab system solver for now)

A2=full(A); %gets spy graph of sparse matrix A
spy(A2)
title("2d Poisson Matrix, n=10")

%plot_calc_solution_P2D(U,grid_x,grid_y,n) %plot solution calculated with discete finite differences
plot_real_solution_P2D(x1,x2,y1,y2,Actual) %plot real soloution 

[max_error,average_error]=error_P2D(Actual,U,x1,x2,y1,y2,grid_x,grid_y,0); %calcs and plots error (last input boolean to determine if plot results)
%% SHows sublinear relationship between error and n
% count=100; 
% max_error=zeros(count,1);
% average_error=zeros(count,1);
% Ns=zeros(count,1);
% for i =5:5+count
%     Ns(i)=i;
%     [A,rhs,grid_x,grid_y]=create_system_P2D(i,x1,x2,y1,y2,G,F); %creates sparse banded matrix to solve problem
%     U=A\rhs; %solves system 
%     [max_error(i),average_error(i)]=error_P2D(Actual,U,x1,x2,y1,y2,grid_x,grid_y,0);
% end
% figure()
% hold on;
% plot(Ns,log10(average_error))
% plot(Ns,log10(max_error))
% title("Max and Average Absolute Error vs n")
% xlabel("n")
% ylabel("Log Error")
% legend("Average Error","Maximum Error")


max_n=10; %max size of n for different solver testing
spacing=2; %step size of n between tests (must be divisor of max_n)
create_system_simp=@(i)create_system_P2D(i,x1,x2,y1,y2,G,F); %function handle to spped up test of same equation at different n values

%method_test(max_n,spacing,create_system_simp) %test differetn system solving methods and plot time vs size of matrix

function []=method_test(max_n,spacing,create_system_simp)
    times=NaN(max_n/spacing,3);
    test_n=zeros(max_n/spacing,1);
    for i=spacing:spacing:max_n %loops through different values for n to compare times
        [A,rhs,~,~]=create_system_simp(i); %create A and rhs to solve
        test_n(i/spacing)=i;
        if i<=80 %stop using regulular matrix in n>80
            tic;
            U=full(A)\rhs; %solves system, default matlab system solver (in documentation looks like this is LU)
            times(i/spacing,1)=toc;
        end
        if i<=24 %stop trying regular solver once n>20
            tic;
            U2=rref(cat(2,full(A),rhs)); %solves system using normal gaussian elimination
            U2=U2(:,i^2+1);
            times(i/spacing,2)=toc;
        end
        tic;
        U_sparse=A\rhs;  %solves it using matlab built in sparse solver (looks like banded solver)
        times(i/spacing,3)=toc;
    end
    %% 
    
    figure() %plots the solve times for eahc method vs n (size n^2 matrix)
    hold on;
    plot(test_n,times(:,1))
    plot(test_n,times(:,2))
    plot(test_n,times(:,3))
    title("Time to Solve Linear System versus n (size n^2 matrix)")
    xlabel("n")
    ylabel("Solve Time (s)")
    legend("Matlab Default (backslash)","Partial Pivoting", "Default Matlab Sparse Solver")
end

function [z]=f(x,y) %rhs forcing function
    z=-8*pi^2*sin(2*pi*y).*sin(2*pi*x);
end

function [z]=g(x,y) %boundry condition function
    z=x.*y; % zero for dirichlet boundry conditions
end

function [z]=actual(x,y) %actual solution to ivp (for plots and error calcs)
    z=sin(2*pi*x).*sin(2*pi*y)+x.*y;
end

function [A,rhs,grid_x,grid_y]=create_system_P2D(n,x1,x2,y1,y2,g,f) %Creates rhs and matrix of a linear system, 
    % also returns grid of poins used in discretization

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
    
    A=spalloc(n^2,n^2,5*n^2); %creates empty variables
    rhs=zeros(n^2,1);
    for i=1:n^2
        if i<=n || i>n^2-n || mod(i,n)==0 || mod(i,n)==1  %checks if u_i is on boundry
            rhs(i)=g(D(i,1),D(i,2)); %if yes, rhs is boundry function
            A(i,i)=1; %corresponding row in A is just row of identity matrix
        else
            rhs(i)=f(D(i,1),D(i,2)); %if no, rhs is f
            A(i,i)=-4/h^2; %approximates second derivvative in both directions at point
            A(i,i+1)=1/h^2;
            A(i,i-1)=1/h^2;
            A(i,i+(n))=1/h^2;
            A(i,i-(n))=1/h^2;
        end
    end
end

function []=plot_calc_solution_P2D(U,grid_x,grid_y,n) %plots calculated solutionq
    grid_u=zeros(n,n); %converts solution back to matrix form
    for i=1:n
        for j=1:n
            grid_u(i,j)=U(n*(i-1)+j);
        end
    end
    figure()
    hold on;
    surface(grid_x,grid_y,grid_u,"FaceColor",'#0072BD',"FaceAlpha",1) %plots solution
end

function []=plot_real_solution_P2D(x1,x2,y1,y2,actual) %plots real solution to BVP
    % calcs and plot actual solution
    [real_x,real_y]=meshgrid(x1:(x2-x1)/150:x2,y1:(x2-x1)/150:y2);
    surface(real_x,real_y,actual(real_x,real_y),"FaceColor",'r','EdgeColor', 'interp',"FaceAlpha",0.3)
    view(3)
end

function [max_err,avg_err]=error_P2D(actual,U,x1,x2,y1,y2,grid_x,grid_y,plot_results) %error cals and plots
    V=scatteredInterpolant(reshape(grid_x,[],1),reshape(grid_y,[],1),U); %creates function that inteprolates calculated solution at V(x,y)
    
    c=200; %number of points to check in interval
    x=linspace(x1,x2,c);
    y=linspace(y1,y2,c);
    solution_error=zeros(c,c);

    for i=1:c
        for j=1:c
            solution_error(i,j)=abs(V(x(i),y(j))-actual(x(i),y(j))); %evalutes relative error of solution at eahc point
        end
    end
    max_err=max(solution_error,[],"all"); %finds max error
    avg_err=mean2(solution_error); %finds average error

    if plot_results==1 %check if result plot is wanted
        [error_x,error_y]=meshgrid(x1:(x2-x1)/(c-1):x2,y1:(x2-x1)/(c-1):y2); %creates grid of x and y value to evaluete error at
        figure()
        surface(error_x,error_y,log10(solution_error)) %plots log10 error
        title("Log of Absolute Error")
        view(3)
    end
end

