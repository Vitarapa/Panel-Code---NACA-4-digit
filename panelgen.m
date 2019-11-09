%This function generates a NACA 4-series airfoil of unit chord length and
%discretizes it into N panels. An additional (N+1)th panel is added to model 
%the airfoil's wake. 
%
%Syntax
%[x,z]=panelgen(airf,N,AoA)
%
%Inputs
%airf - 4-character string  - used to decribe the airfoil characteristics.
%                             For example if airf = '2412' a NACA 2412
%                             airfoil is generated
%N    - Integer scalar      - Number of panels to be generated
%AoA  - Real scalar         - Angle of attack of airfoil (alpha) in radians
%
%Outputs
%x, z - real 1D arrays of length N+2 - the arrays contain the coordinates 
%                                      of all panel endpoints. The endpoints
%                                      for the ith panel are stored in the 
%                                      ith and (i+1)th cells 

function [x,z]=panelgen(airf,N,AoA)

%retreive information from NACA code
if length(airf)==4
    c=str2num(airf(1))/100; %read max camber
    p=str2num(airf(2))/10; %read position of max camber
    t=str2num(airf(3:4))/100; %read thickness to chord ratio
else
    error('The airfoil descriptor should have 4 characters')
end
   
x=zeros(N+2,1);
z=x;
for i=1:N+1
    xt=1-0.5*(1-cos(2*pi*(i-1)/N)); %calculate panel endpoint x using cos distribution
    %calculate y ordinate for symmetric airfoil
    if i==1 || i==N+1
        yt=0; %ensure airfoil comes back to a point.
    else
        yt=5*t*(0.2969*sqrt(xt)-0.126*xt-0.3516*xt^2+0.2843*xt^3-0.1015*xt^4);
    end
    %calculate cambeline height and angle
    if (c>0)
        if xt<=p %if forward of max caber point
            theta=c*2*(p-xt)/p^2;
            yc=c*(2*p*xt-xt^2)/p^2;
        else %if aft of max camber point
            theta=c*2*(p-xt)/(1-p)^2;
            yc=c*((1-2*p)+2*p*xt-xt^2)/(1-p)^2;
        end
    else %if no camber
        theta=0;
        yc=0;
    end
    %apply camber and get final airfoil shape
    if i-1>=0.5*N
        x(i)=xt-yt*sin(theta);
        z(i)=yc+yt*cos(theta);
    else
        x(i)=xt+yt*sin(theta);
        z(i)=yc-yt*cos(theta);
    end
end
%add wake panel endpoint
x(N+2)=10000*cos(AoA);
z(N+2)=10000*sin(AoA);
    