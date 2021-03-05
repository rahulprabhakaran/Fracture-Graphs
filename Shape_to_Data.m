function [ apodi_2, xmin, xmax, ymin, ymax ] = Shape_to_Data( B )
% converts fracture shape file into an array of edges
% R Prabhakaran (2017)
[m1,~]=size(B);
for i=1:m1
     [n,~]=size((B(i).X)');
     frac_index(i,1)=n-1;
end
frac_index(:,2)=cumsum(frac_index(:,1));
frac_index(:,3)=frac_index(:,1)-1;
frac_index(:,4)=cumsum(frac_index(:,3));

M=max(frac_index(:,4));
apodi_2=zeros(M,6); % change 5 -> 6 if there is more data
x=1;

for i=1:m1
    m2=frac_index(i,1);
    m3=frac_index(i,2);
    m4=frac_index(i,3);
      for j=1:m2-1
        apodi_2(x,1)=B(i).X(j);
        apodi_2(x,2)=B(i).Y(j); 
        apodi_2(x,3)=B(i).X(j+1);
        apodi_2(x,4)=B(i).Y(j+1);
        %apodi_2(x,5)=B(i).Layer;
        %apodi_2(x,5)=B(i).Polyline_1;
        %apodi_2(x,5)=B(i).Polyline_ID;
%           apodi_2(x,5)=B(i).Tile_ID;
%            apodi_2(x,5)=B(i).fracture_i;
%        apodi_2(x,5)=B(i).FracID;
        %apodi_2(x,5)=B(i).Id;
%        apodi_2(x,6)=B(i).LENGTH;
        x=x+1;
      end     
end

frac_id=unique(apodi_2(:,5),'rows');
N=length(frac_id);
for i=1:N
  leng_B(i,1) = nnz(apodi_2(:,5)==frac_id(i,1));  
end
leng_B(:,2)=cumsum(leng_B(:,1));
N1N2_Y=zeros(N,2);
N1N2_Y(1,1)=1;
N1N2_Y(1,2)=12;
for i=2:N
 N1N2_Y(i,1)=leng_B(i-1,2)+1; 
 N1N2_Y(i,2)=leng_B(i,2);
 if N1N2_Y(i,1)>N1N2_Y(i,2)
    N1N2_Y(i,1)= N1N2_Y(i,2);
 end
end
N1N2_Y(:,3)=apodi_2(N1N2_Y(:,2),6);
for i=1:N
 Data_Raw(i,1:2)=apodi_2(N1N2_Y(i,1),1:2);
 Data_Raw(i,3:4)=apodi_2(N1N2_Y(i,2),3:4);
end

xmin = min([apodi_2(:,1);apodi_2(:,3)] );
xmax = max([apodi_2(:,1);apodi_2(:,3)] );
ymin = min([apodi_2(:,2);apodi_2(:,4)] );
ymax = max([apodi_2(:,2);apodi_2(:,4)] );

% for i=1:length(N1N2_Y) 
%  fraclength = sum(Lengths2D(apodi_2(N1N2_Y(i,1):N1N2_Y(i,2),:)));
%  apodi_2(N1N2_Y(i,1):N1N2_Y(i,2),6) = fraclength;
%  clearvars fraclength
% end


end

