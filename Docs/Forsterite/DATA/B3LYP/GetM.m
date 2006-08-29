
function [M]=GetM(t)
# t=[a b c alpha beta gamma]
#Compute M=(a,b,c)
DegToRad=pi/180;
M=zeros(3,3);
M(1,1)=t(1);
M(1,2)=t(2)*cos(DegToRad*t(6));
M(2,2)=t(2)*sin(DegToRad*t(6));
M(1,3)=t(3)*cos(DegToRad*t(5));
M(2,3)=(t(2)*t(3)*cos(DegToRad*t(4))-M(1,2)*M(1,3))/M(2,2);
M(3,3)=sqrt(t(3)^2-M(1,3)^2-M(2,3)^2);


