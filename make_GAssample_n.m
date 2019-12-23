function [ TW ] = make_GAssample_n( parameter )

ge=1;

sam=1;
% parameter=lhsdesign(sam,3*ge);
%
% for ii=1:3*ge
%     if rem(ii,3)==1
%         parameter(:,ii)=6*parameter(:,ii);
%     elseif rem(ii,3)==2
%         parameter(:,ii)=6*parameter(:,ii);
%     elseif ii==3;
%         parameter(:,ii)=0.1*parameter(:,ii);
%     else
%         parameter(:,ii)=2*parameter(:,ii)-1;
%     end
% end
sp1=100; % partitioned element number
sp2=50; % partitioned element number for the opposite side
ne=2*sp1*sp2+sp2^2; % number of elements

for ll=1:sam
    elementcordinate=load('C:\Users\dora-user\Documents\MATLAB\elm_hashin.txt');
    
    X=elementcordinate(:,1)-0.2;
    Y=elementcordinate(:,2)-0.4;
    
    temp = Y;
    Y = -X;
    X = temp;
    
    pp=1:3*ge;
    for j=1:3*ge
        pp(j)=parameter(ll,j);
    end
    
    xb=1:4*ge;
    yb=1:4*ge;
    for j=1:ge
        xb(4*j-3)=x(4*j-3)*a^2/(x(4*j-3)^2+y(4*j-3)^2);
        yb(4*j-3)=y(4*j-3)*a^2/(x(4*j-3)^2+y(4*j-3)^2);
        xb(4*j-2)=x(4*j-2)*a^2/(x(4*j-2)^2+y(4*j-2)^2);
        yb(4*j-2)=y(4*j-2)*a^2/(x(4*j-2)^2+y(4*j-2)^2);
        xb(4*j-1)=x(4*j-1)*a^2/(x(4*j-1)^2+y(4*j-1)^2);
        yb(4*j-1)=y(4*j-1)*a^2/(x(4*j-1)^2+y(4*j-1)^2);
        xb(4*j)=x(4*j)*a^2/(x(4*j)^2+y(4*j)^2);
        yb(4*j)=y(4*j)*a^2/(x(4*j)^2+y(4*j)^2);
    end 
    
    U=zeros(ne,1);
    V=zeros(ne,1);
    for j=1:ge
        U = U -K(j)*(X-x(4*j-3))./((X-x(4*j-3)).^2+(Y-y(4*j-3)).^2) -K(j)*(X-x(4*j-2))./((X-x(4*j-2)).^2+(Y-y(4*j-2)).^2) +K(j)*(X-x(4*j-1))./((X-x(4*j-1)).^2+(Y-y(4*j-1)).^2) +K(j)*(X-x(4*j))./((X-x(4*j)).^2+(Y-y(4*j)).^2)...
            -K(j)*(X-xb(4*j-3))./((X-xb(4*j-3)).^2+(Y-yb(4*j-3)).^2) -K(j)*(X-xb(4*j-2))./((X-xb(4*j-2)).^2+(Y-yb(4*j-2)).^2) +K(j)*(X-xb(4*j-1))./((X-xb(4*j-1)).^2+(Y-yb(4*j-1)).^2) +K(j)*(X-xb(4*j))./((X-xb(4*j)).^2+(Y-yb(4*j)).^2)...
            ;
        
        V = V -K(j)*(Y-y(4*j-3))./((X-x(4*j-3)).^2+(Y-y(4*j-3)).^2) -K(j)*(Y-y(4*j-2))./((X-x(4*j-2)).^2+(Y-y(4*j-2)).^2) +K(j)*(Y-y(4*j-1))./((X-x(4*j-1)).^2+(Y-y(4*j-1)).^2) +K(j)*(Y-y(4*j))./((X-x(4*j)).^2+(Y-y(4*j)).^2) ...
            -K(j)*(Y-yb(4*j-3))./((X-xb(4*j-3)).^2+(Y-yb(4*j-3)).^2) -K(j)*(Y-yb(4*j-2))./((X-xb(4*j-2)).^2+(Y-yb(4*j-2)).^2) +K(j)*(Y-yb(4*j-1))./((X-xb(4*j-1)).^2+(Y-yb(4*j-1)).^2) +K(j)*(Y-yb(4*j))./((X-xb(4*j)).^2+(Y-yb(4*j)).^2)...
            ;
    end
    
    fiberangle = atand(V./U) + 90;
    Si=1:ne;
    roop=ne;
    fid = fopen( 'elmangle_hashin.inp', 'w+');
    
    for i = 1:roop
        fprintf(fid, 'S%d=%f\n',Si(i),fiberangle(i,:));
    end
    
    fclose(fid);
    
    fid = fopen( 'inputdata_h_n3.txt', 'a+');
    for j=1:3*ge
        if rem(j,3*ge)~=0
            fprintf(fid,'%f ',pp(j));
        else
            fprintf(fid,'%f\n',pp(j));
        end
    end
    
    fclose(fid);
    system('SET KMP_STACKSIZE=2048k & "E:\Program Files\ANSYS Inc\v171\ansys\bin\winx64\ANSYS171.exe" -b aa_t_a -i "E:\ansys\experimental_p2.txt" -o "E:\ansys\out.txt"');
    TW = load('C:\Users\dora-user\Documents\MATLAB\Hashin.txt');
    fid = fopen( 'outputdata_h_n3.txt', 'a+');
    fprintf(fid,'%f\n',TW);
    fclose(fid);
   
end
end
