%Function to compute the u and v velocity components produced at a point p with
%coordinates p=[xp,zp] by a constant strength doublet panel with end points 
%at p1=[x1,z1] and p2=[x2,z2] and unit strenth.
%
%the function syntax is [u,v]=cdoublet(p,p1,p2)
%
%Outputs
%u - real number - horizontal velocity component due to doublet in global FoR
%v - real number - vertical velocuty component due to doublet in global FoR
%
%Inputs
%p  - 1 x 2 array - coordinates [x,z] of point at which to evaluate u,v
%p1 - 1 x 2 array - coordinates [x1,z1] of panel start point
%p2 - 1 x 2 array - coordinates [x2,z2] of panel end point
%mu - real number - doublet strength

function [u,v]=cdoublet2(p,p1,p2)

%convert coordinate system to doublet frame of reference
alfa=-atan2(p2(:,2)-p1(:,2),p2(:,1)-p1(:,1)); %find local panel angle
T=[cos(alfa(1)),-sin(alfa(1));sin(alfa(1)),cos(alfa(1))];%transofrmation matrix
A=T*(p-p1)';
B=T*(p-p2)';
d1=A(1,:); %x-x1 in local FoR
d2=B(1,:); %x-x2 in local FoR
dz=A(2,:); %height above doublet in local FoR


u1 = zeros(1,length(dz));
v1 = zeros(1,length(dz));
%estimate u,v in doublet FoR
index = find(abs(dz)<0.000001);
% u1(index) = 0;
v1(index) = (d1(index)./((d1(index)).^2)-d2(index)./(d2(index).^2))/(2*pi);

index2 = find(abs(dz)>=0.000001);
u1(index2)=-dz(index2).*(1./(d1(index2).^2+dz(index2).^2)-1./(d2(index2).^2+dz(index2).^2))./(2*pi);
v1(index2)=(d1(index2)./(d1(index2).^2+dz(index2).^2)-d2(index2)./(d2(index2).^2+dz(index2).^2))/(2*pi);
% if abs(dz)<0.000001
%     u1=0;
%     v1=(d1/(d1^2)-d2/(d2^2))/(2*pi);
% else
%     u1=-dz*(1/(d1^2+dz^2)-1/(d2^2+dz^2))/(2*pi);
%     v1=(d1/(d1^2+dz^2)-d2/(d2^2+dz^2))/(2*pi);
% end
%     
%convert u,v to gloval Frame of Reference
temp=[cos(alfa(1)),sin(alfa(1));-sin(alfa(1)),cos(alfa(1))]*[u1;v1];
u=temp(1,:)';
v=temp(2,:)';


