function lap=laplacian2(dx,dy,u)
    n1=length(u(:,1));
    n2=length(u(1,:));
    lap=zeros(n1,n2);
    for i=1:n1
        for j=2:n2-1
           if i==1
               lap(i,j)=(-2*u(i,j)+2*u(i+1,j))/dx^2+(u(i,j-1)-2*u(i,j)+u(i,j+1))/dy^2;
           elseif i==n1
               lap(i,j)=(-2*u(i,j)+2*u(i-1,j))/dx^2+(u(i,j-1)-2*u(i,j)+u(i,j+1))/dy^2;
           else 
           lap(i,j)=(u(i-1,j)-2*u(i,j)+u(i+1,j))/dx^2+(u(i,j-1)-2*u(i,j)+u(i,j+1))/dy^2;
           end
        end
    end
   
   for i=1:n1
           if i==1
               lap(i,1)=(-2*u(i,1)+2*u(i+1,1))/dx^2+(-2*u(i,1)+2*u(i,2))/dy^2;
               lap(i,n2)=(-2*u(i,n2)+2*u(i+1,n2))/dx^2+(-2*u(i,n2)+2*u(i,n2-1))/dy^2;
           elseif i==n1
               lap(i,1)=(-2*u(i,1)+2*u(i-1,1))/dx^2+(-2*u(i,1)+2*u(i,2))/dy^2;
               lap(i,n2)=(-2*u(i,n2)+2*u(i-1,n2))/dx^2+(-2*u(i,n2)+2*u(i,n2-1))/dy^2;
           else 
           lap(i,1)=(u(i-1,1)-2*u(i,1)+u(i+1,1))/dx^2+(-2*u(i,1)+2*u(i,2))/dy^2;
           lap(i,n2)=(u(i-1,n2)-2*u(i,n2)+u(i+1,n2))/dx^2+(-2*u(i,n2)+2*u(i,n2-1))/dy^2;
           end
   end
end