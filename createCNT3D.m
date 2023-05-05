%%This code is developed to generate CNT/graphene junction.
%%In this code, carbon nanotube is allowed to rotate along 3 directions.
%%Developed in 2017
tic;
clear;clc;
close all;
dbstop if error
global n m a;%%n and m are carbon nanotube indices, a is bond length.
n=5;  %%n>==m
m=0;
switch_guess_plot=0; %%plot guess(0/1);
switch_structure_plot=0;%%plot pre-structure(0/1);
switch_reslt_plot=1; %%plot result(0/1);
switch_area_select=0; %%graphene area:0--circle,1--ring;
switch_polygon=0; %% non-hexagon defects except heptagon and pentagon: forbid--1,allow--0 Only count the number of polygons, need further development to control structure.
z=11;%%CNT length
D_cutoff=1.95; %%unit in angstrom
obj=['CNT(',num2str(n),',',num2str(m),',',num2str(z),').dat'];%%generate graphene data
delete(obj);

a=1.42; %%unit in angstrom
if m==n
    zz=2*a*z*cos(pi/6);
elseif mod(z,2)==1
    zz=2*a+(z-1)/2*3*a;
else
    zz=a*(sin(pi/6)+z*3/2);
end

%decide the rectangle to form CNT
[CD,x1,x2,y1,y2]=Rec(m,n,zz); %decide the rectangle to form CNT
L = ceil(max([abs(CD(1,1)), abs(x1),abs(x2)])/(2*a*cos(pi/6)))+1;
HH= ceil(max([abs(CD(1,2)), abs(y1),abs(y2)])/(2.5*a));
H=4*HH;

y(1,1) = 0;
x(1,1) = a*cos(pi/6);
flag = 0;
H=70;
L=30;

for i = 2:H
    if mod(i,2)==1
        y(1,i) = y(1,i-1) + a;
    else
        y(1,i) = y(1,i-1) + a/2;
        x(1,i) = flag;
        if (i+1)<= H
            x(1,i+1) = flag;
        end
        if flag == 0
            flag = a*cos(pi/6);
        else
            flag = 0;
        end
    end
end


for i = 2: L
    for j = 1: H
        x(i,j) = x(i-1,j) + 2*a*cos(pi/6);
        y(i,j) = y(1,j);
    end
end

x = x-a*cos(pi/6); 
figure (1)

xf = x;
yf = -y-a;
xt=-xf;
yt=yf;

xx = [x,xf];
yy = [y,yf];
[q,p] = size(x);


for i = 1:L
    plot (x(i,:),y(i,:), 'ko','MarkerFaceColor','k');
    line(x(i,:),y(i,:),'color','k','linewidth',4)
    for j= 1:H
        if i<L&&j<=floor(p/4)
            line([x(i,4+(j-1)*4),x(i+1,3+(j-1)*4)],[y(i,4+(j-1)*4),y(i+1,3+(j-1)*4)],'color','k','linewidth',4)
            line([x(i,1+(j-1)*4),x(i+1,2+(j-1)*4)],[y(i,1+(j-1)*4),y(i+1,2+(j-1)*4)],'color','k','linewidth',4)
        end
    end
    hold on;
end
axis equal


recxv=[x(1,1), CD(1, 1), x1+1e4*eps, x2, x(1,1)];
recyv=[y(1,1), CD(1, 2)+1e4*eps, y1+1e4*eps, y2-1e4*eps, y(1,1)];
plot(recxv,recyv,'b','LineWidth',2); 
[in,on] = inpolygon(xx,yy,recxv,recyv);
plot(xx(in),yy(in),'go')
oy=[0,1];
[xnew,ynew]=tl(xx(in),yy(in),CD,oy);%% rotate the retangle and make a side on the X axis
r=norm(CD)/(2*pi);%% carbon nanotube radius

%%roll the graphene to form carbon nanotube
ynew_max=max(ynew);
l=abs(ynew-ynew_max)<1e2*eps;
xnew(l)=[];
ynew(l)=[];
theta=ynew/r;

figure (2)%% plot the carbon nanotube bottom in polar
polar(theta,r*ones(size(theta)),'x')
[xxx,yyy,zzz] = pol2cart(theta,r*ones(size(theta)),xnew);
b=[xxx;yyy;zzz]';
b(abs(b)<1e-13)=0;
c=unique(b,'rows','stable');
c=c';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Generate two data strucutres%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (m~=n) && (m~=0)
    n_bottom_c=length(find(c(3,:)<4*a));
    any_two_distance_c=zeros(n_bottom_c);
    num_bottom_c=find(c(3,:)<4*a);
    
    neighbor_list_c=zeros(n_bottom_c);
    %
    for nx = 1:n_bottom_c
        for ny = 1:n_bottom_c
            
            if nx==ny
                continue;
            end
            
            temp = norm(c(:, num_bottom_c(nx)) - c(:, num_bottom_c(ny)));
            any_two_distance_c(nx,ny)=temp;
            
            
        end
    end
    value_any_two_distance_c=unique(any_two_distance_c);
    for set_value=1:length(value_any_two_distance_c)
        if value_any_two_distance_c(set_value)~=0
            value_any_two_distance_c_1=value_any_two_distance_c(set_value);
            break
        end
        
    end
    for set_value=1:length(value_any_two_distance_c)
        if value_any_two_distance_c(set_value)>value_any_two_distance_c_1+0.3
            value_any_two_distance_c_2=value_any_two_distance_c(set_value);
            break
        end
    end
    for set_value=1:length(value_any_two_distance_c)
        if value_any_two_distance_c(set_value)>value_any_two_distance_c_2+0.3
            value_any_two_distance_c_3=value_any_two_distance_c(set_value);
            break
        end
    end
    
    for nx = 1:n_bottom_c
        for ny = 1:n_bottom_c
            
            if nx==ny
                continue;
            end
            
            if (any_two_distance_c(nx,ny) > value_any_two_distance_c_1-0.1) && (any_two_distance_c(nx,ny) < value_any_two_distance_c_1+0.3)
                neighbor_list_c(nx,ny)=1;
            elseif (any_two_distance_c(nx,ny) > value_any_two_distance_c_2-0.1) && (any_two_distance_c(nx,ny) < value_any_two_distance_c_2+0.3)
                neighbor_list_c(nx,ny)=2;
            elseif (any_two_distance_c(nx,ny) > value_any_two_distance_c_3-0.1) && (any_two_distance_c(nx,ny) < value_any_two_distance_c_3+0.3)
                neighbor_list_c(nx,ny)=3;
            end
            
        end
    end
    del_atom_cnt=zeros(1,n_bottom_c);
    bottom_atom_cnt=zeros(1,n_bottom_c);

    for ny = 1:n_bottom_c
        if c(3,num_bottom_c(ny))<a
            if length(find(neighbor_list_c(:,ny)==1))==1;
                del_atom_cnt(ny)=1;
            end
        end
        if c(3,num_bottom_c(ny))<2*a
            if length(find(neighbor_list_c(:,ny)==1))==2;
                bottom_atom_cnt(ny)=1;
            end
        end
    end

    neighbor_list_c(:,find(del_atom_cnt==1))=0;
    neighbor_list_c(find(del_atom_cnt==1),:)=0;
    for ny = 1:n_bottom_c
        if c(3,num_bottom_c(ny))<a
            if length(find(neighbor_list_c(:,ny)==1))==1;
                del_atom_cnt(ny)=1;
            end
        end
        if c(3,num_bottom_c(ny))<2*a
            if length(find(neighbor_list_c(:,ny)==1))==2;
                bottom_atom_cnt(ny)=1;
            end
        end
    end
    
    c_bottom_atom_cnt=c(:,num_bottom_c(find(bottom_atom_cnt==1)));
    c(:,num_bottom_c(find(del_atom_cnt==1)))=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xxx=c(1,:);
yyy=c(2,:);
zzz=c(3,:);


figure (4)%% plot the carbon nanotube
plot3(xxx,yyy,zzz,'o')
for i = 1:length(c)
    text(c(1,i)+0.05, c(2,i)+0.05,c(3,i)+0.05, int2str(i), 'tag', 'text');
end
axis equal

%% output carbon nanotube file for ovito
filename=['CNT(',num2str(n),',',num2str(m),',',num2str(z),').dat'];
fid=fopen(filename,'a');
leng=length(xxx);
Out1=[num2str(leng),' atoms\n'];
Out2=[num2str(min(xxx)-1.5),' ',num2str(max(xxx)+1.5),' xlo xhi\n'];
Out3=[num2str(min(yyy)-1.5),' ',num2str(max(yyy)+1.5),' ylo yhi\n'];
Out4=[num2str(min(zzz)-1.5),' ',num2str(max(zzz)+1.5),' zlo zhi\n'];
fprintf(fid,'# LAMMPS data file\n');
fprintf(fid,Out1);
fprintf(fid,'1 atom types\n');
fprintf(fid,Out2);
fprintf(fid,Out3);
fprintf(fid,Out4);
fprintf(fid,'\n');
fprintf(fid,'Atoms\n');
fprintf(fid,'\n');
atom=1;
for j=1:1:leng
    fprintf(fid,'%d %d %12.6f %12.6f %12.6f\n', j, atom, xxx(j), yyy(j), zzz(j));
end

fclose(fid);


%Carbon nanotube bottom
if (m~=n & m~=0)
    cc=c_bottom_atom_cnt;
else
    cc=c;
    ll=(cc(3,:)>1e3*eps);
    cc(:,ll)=[];
end
figure (5)
plot(cc(1,:),cc(2,:),'rx')
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select part of the graphene to reduce calculation

H_prime=70;
L_prime=30;
y_prime(1,1) = 0;
x_prime(1,1) = a*cos(pi/6);
flag = 0;


for i = 2:H_prime
    if mod(i,2)==1
        y_prime(1,i) = y_prime(1,i-1) + a;
    else
        y_prime(1,i) = y_prime(1,i-1) + a/2;
        x_prime(1,i) = flag;
        if (i+1)<= H_prime
            x_prime(1,i+1) = flag;
        end
        if flag == 0
            flag = a*cos(pi/6);
        else
            flag = 0;
        end
    end
end


for i = 2: L_prime
    for j = 1: H_prime
        x_prime(i,j) = x_prime(i-1,j) + 2*a*cos(pi/6);
        y_prime(i,j) = y_prime(1,j);
    end
end

x_prime = x_prime-a*cos(pi/6);
xf_prime = x_prime;
yf_prime = -y_prime-a;
xt=-xf_prime;
yt=yf_prime;

xx_prime = [x_prime,xf_prime];
yy_prime = [y_prime,yf_prime];


xxx_prime=xx_prime(:);
yyy_prime=yy_prime(:);
graphene_prime=unique([xxx_prime, yyy_prime],'rows');
graphene_prime=[graphene_prime(:,1)-2.459512146747806,graphene_prime(:,2)];
graphene_prime(find(graphene_prime(:,1)<-1 | graphene_prime(:,1)>50),:)=[];
graphene_prime(find(graphene_prime(:,2)<-26 | graphene_prime(:,2)>27),:)=[];
graphene_prime=[graphene_prime(:,1)-24.595121467478055,graphene_prime(:,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sgraph=graphene_prime'; %% coordinates of selected graphene atoms


sltr_zigzap=sqrt((cc(1,1)-cc(2,1))^2+(cc(2,2)-cc(2,2))^2)/2;
sizecc=length(cc);

%%Set atom ratio that allow the graphene atoms inside the circle composed
%%of the cnt bottom atoms
if sizecc>=8
    tolerance=0.9;
else
    tolerance=0.7;
end

%%Set atom ratio that make sure the distance between the graphene atom and
%%paired CNT atom is smaller than some certain value
if sizecc>=10
    sizecc_partial=sizecc*2/3;
else
    sizecc_partial=sizecc/2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select atoms on graphene to bond with carbon nanotube%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N_neighbor;
global neighbor_list scgraph slct_atom_index any_two_distance ssgraphatom_0;
%%neighbor type for any two graphene atoms, coordinates of graphene atoms, 
%%index of selected graphene atoms, distance between any two graphene atoms
%%loop starting point of selected graphene atoms

if switch_area_select==0
    scgraph_index=find(((sgraph(1,:).^2+(sgraph(2,:).^2)<(norm(cc(1:2,1))+5*1.42)^2)));
else
    f_res_area=sqrt(n^2+m^2+m*n);
    d_res_area=sqrt(3)*a/pi*f_res_area;
    if f_res_area>=6
        scgraph_index=find((sgraph(1,:).^2+(sgraph(2,:).^2)<(d_res_area/2+3*a)^2) & (sgraph(1,:).^2+(sgraph(2,:).^2)>(d_res_area/2-3*a)^2));
    else
        scgraph_index=find((sgraph(1,:).^2+(sgraph(2,:).^2)<(d_res_area/2+3*a)^2));      
    end
end

scgraph=sgraph(:,scgraph_index);

figure (8) %% plot the graphene atoms (black) and carbon nanotue bottom atoms(red)
plot(scgraph(1,:),scgraph(2,:),'k.','MarkerSize',30)
hold on
plot(cc(1,:),cc(2,:),'r.','MarkerSize',30)
axis equal

distance_0=sqrt((scgraph(1,:)-cc(1,1)).^2 + (scgraph(2,:)-cc(2,1)).^2);
ssgraphatom_ind_0=find(distance_0==min(distance_0)); % initial point
ssgraphatom_0=scgraph(:,ssgraphatom_ind_0(1));
scgraph_number=length(scgraph_index);
plot(ssgraphatom_0(1,:),ssgraphatom_0(2,:),'g.','MarkerEdgeColor','g','MarkerSize',30)
for i = 1:scgraph_number
    text(scgraph(1,i)+0.25, scgraph(2,i)+0.05, int2str(i), 'tag', 'text');
end

hold off


any_two_distance=zeros(scgraph_number);
neighbor_list=zeros(scgraph_number);

%% obtain neighbor type between any two graphene atoms
for x = 1:scgraph_number
    for y = 1:scgraph_number
        
        if x==y
            continue;
        end
        
        temp = norm(scgraph(:, x) - scgraph(:, y));
        any_two_distance(x,y)=temp;
        if temp<a+eps*1e2
            neighbor_list(x,y)=1;
        elseif temp<sqrt(3)*a+eps*1e2
            neighbor_list(x,y)=2;
        elseif temp < 2*a+eps*1e2
            neighbor_list(x,y)=3;
        end
    end
end
N_neighbor=scgraph_number-sum(~neighbor_list,2);

slct_atom_index=zeros(1,n+m);
slct_atom_index(1)=ssgraphatom_ind_0;

current_num=1;
line=1;
output_N=2000;
output=zeros(output_N,n+m);
[output,line]=guess_next(ssgraphatom_ind_0, current_num,line,output);

if line<output_N
    output(line:output_N,:)=[];
end
if switch_guess_plot==1
    for image_number=1:length(output)
        figure(100+image_number)%%plot all possible group of selected graphene atoms
        plot(scgraph(1,:),scgraph(2,:),'*');
        hold on
        plot(scgraph(1,output(image_number,:)),scgraph(2,output(image_number,:)),'bo');
        hold off
        pause;
    end
end
output_prime=[output,output(:,1)];
output_copy=output;
output_neigh_num=zeros(size(output));
[length_output,width_output]=size(output);
for i_check_slct=1:length_output
    for j_check_slct=1:width_output
        output_neigh_num(i_check_slct,j_check_slct)=neighbor_list(output_prime(i_check_slct,j_check_slct+1),output_prime(i_check_slct,j_check_slct));
    end
end
output_neigh_num_prime=output_neigh_num;
output_neigh_num_prime(output_neigh_num_prime(:)~=1)=0;
temp_check_slct1=zeros(1,n+m);
temp_check_slct_index=zeros(1,length_output);
for i_onnp=1:length_output
    if sum(output_neigh_num_prime(i_onnp,:))<=1
        continue
    end
    temp_check_slct1=find(output_neigh_num_prime(i_onnp,:)==1);
    diff(temp_check_slct1);
    
    if (sum(find(diff(temp_check_slct1)==1))~=0)
        temp_check_slct_index(i_onnp)=1;
        continue
    end
    if (output_neigh_num_prime(i_onnp,1)==1 && output_neigh_num_prime(i_onnp,n+m)==1)
        temp_check_slct_index(i_onnp)=1;
        continue
    end
    
end
output_neigh_num(find(temp_check_slct_index==1),:) = [];
output_copy(find(temp_check_slct_index==1),:) = [];

coordin_output_copy=zeros(2,n+m);
[length_output_copy,width_output_copy]=size(output_copy);
temp_check_sp2_index=zeros(1,length_output_copy);
for i_onnp2=1:length_output_copy
    coordin_output_copy=scgraph(:,output_copy(i_onnp2,:));
    scgraph_copy=scgraph;
    [in_scgraph_copy,on_scgraph_copy] = inpolygon(scgraph_copy(1,:),scgraph_copy(2,:),coordin_output_copy(1,:),coordin_output_copy(2,:));
    del_scgraph_copy=in_scgraph_copy~=on_scgraph_copy;
    neighbor_list_copy=neighbor_list;
    neighbor_list_copy(del_scgraph_copy,:)=0;
    neighbor_list_copy(:,del_scgraph_copy)=0;
    n_sp2_check=n+m;
    for sp2_check=1:n_sp2_check
        sp2_index=find(neighbor_list_copy(output_copy(i_onnp2,sp2_check),:)==1);
        if length(sp2_index)~=2
            temp_check_sp2_index(i_onnp2)=1;
            continue
        end
    end
end
output_neigh_num(find(temp_check_sp2_index==1),:) = [];
output_copy(find(temp_check_sp2_index==1),:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% REPEAT  CHECK%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[i_num,j_CNT]=size(output_neigh_num);
temp_rpt_array=zeros(1,j_CNT);

index_dlt_array=zeros(1,i_num);
for i_rpt_check=1:i_num
    index_rpt_array=[];
    if index_dlt_array(i_rpt_check)==1
        continue
    end
    temp_rpt_array=output_neigh_num(i_rpt_check,:);
    index_rpt_array=[index_rpt_array,find(ismember(output_neigh_num, temp_rpt_array,'rows')==1)];
    
    index_dlt_array(index_rpt_array)=1;
    if length(temp_rpt_array)-length(unique(temp_rpt_array))==j_CNT-1
        index_dlt_array(i_rpt_check)=0;
        continue
    end
    for j_rpt_check=2:j_CNT
        test_temp_rpt_array=[temp_rpt_array(j_rpt_check:j_CNT),temp_rpt_array(1:j_rpt_check-1)];
        index_rpt_array=find(ismember(output_neigh_num,test_temp_rpt_array,'rows')==1);
        index_dlt_array(index_rpt_array)=1;
        fliplr_temp_rpt_array=fliplr(test_temp_rpt_array);
        index_rpt_array=find(ismember(output_neigh_num, fliplr_temp_rpt_array,'rows')==1);
        index_dlt_array(index_rpt_array)=1;
    end
    index_dlt_array(i_rpt_check)=0;
    
end
output_neigh_num(find(index_dlt_array==1),:) = [];
output_copy(find(index_dlt_array==1),:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrangement_ssgraph=output_copy;
length_ssgraphcrd=size(arrangement_ssgraph,1);
ssgraph=zeros(2,sizecc);
Xarray=zeros(1,length_ssgraphcrd);
Yarray=zeros(1,length_ssgraphcrd);
alphaarray=zeros(1,length_ssgraphcrd);
betaarray=zeros(1,length_ssgraphcrd);
gammaarray=zeros(1,length_ssgraphcrd);
harray=zeros(1,length_ssgraphcrd);
subfarray=zeros(1,length_ssgraphcrd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%check polygons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Need further development%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_polygon==1
    
    dlt_ssgraph=zeros(length_ssgraphcrd,1);

    
    for  exh=1:length_ssgraphcrd
        neighbor_list_copy_temp=neighbor_list;
        neighbor_list_order=zeros(1,length(arrangement_ssgraph(exh,:)));
        %%%%%%%%%%%%%%%
        ssgraph_temp(1,:)=scgraph(1,arrangement_ssgraph(exh,:));
        ssgraph_temp(2,:)=scgraph(2,arrangement_ssgraph(exh,:));
        scgraph_temp=[scgraph(1,:);scgraph(2,:)];
        [in_cnt_temp0,on_cnt_temp0] = inpolygon(scgraph(1,:),scgraph(2,:),ssgraph_temp(1,:),ssgraph_temp(2,:));
        del_cnt_temp0=in_cnt_temp0~=on_cnt_temp0;
        scgraph_temp(:,del_cnt_temp0)=[];
        find(del_cnt_temp0==1);
        neighbor_list_copy_temp(del_cnt_temp0,:)=0;
        neighbor_list_copy_temp(:,del_cnt_temp0)=0;
        %%%%%%%%%%%%%%%
        for exh_a=1:length(arrangement_ssgraph(exh,:))
            arrangement_ssgraph_plus=[arrangement_ssgraph(exh,:);arrangement_ssgraph(exh,2:end),arrangement_ssgraph(exh,1)];
            neighbor_list_order(exh_a)=neighbor_list_copy_temp(arrangement_ssgraph_plus(1,exh_a),arrangement_ssgraph_plus(2,exh_a));
            if neighbor_list_order(exh_a)==2
                find(neighbor_list_copy_temp(arrangement_ssgraph_plus(1,exh_a),:)==1)
                find(neighbor_list_copy_temp(arrangement_ssgraph_plus(2,exh_a),:)==1)
                N_neighbor_4=sum(ismember(find(neighbor_list_copy_temp(arrangement_ssgraph_plus(1,exh_a),:)==1),find(neighbor_list_copy_temp(arrangement_ssgraph_plus(2,exh_a),:)==1)));
                if N_neighbor_4==0
                    neighbor_list_order(exh_a)=4;
                end
                
            end
            
        end
        if m==0 && n~=m
            if sum(find(neighbor_list_order)==1)~=0
                dlt_ssgraph(exh)=1;
            end
        else
            length_cc=length(cc);
            distance_cc=zeros(1,length_cc);
            cc_neighbor_list=zeros(1,length_cc);
            for exh_b=1:length_ssgraph-1
                distance_cc(exh_b)=sqrt((cc(1,exh_b+1)-cc(1,exh_b))^2+(cc(2,exh_b+1)-cc(2,exh_b))^2+(cc(3,exh_b+1)-cc(3,exh_b))^2);
                if distance_cc(exh_b)<a+eps*1e2
                    cc_neighbor_list(exh_b)=1;
                elseif distance_cc(exh_b)<sqrt(3)*a+eps*1e2
                    cc_neighbor_list(exh_b)=2;
                elseif distance_cc(exh_b) < 2*a+eps*1e2
                    cc_neighbor_list(exh_b)=3;
                end
            end
            distance_cc(length_cc)=sqrt((cc(1,length_cc)-cc(1,1))^2+(cc(2,length_cc)-cc(2,1))^2+(cc(3,length_cc)-cc(3,1))^2);
            if distance_cc(length_cc)<a+eps*1e2
                cc_neighbor_list(exh_b)=1;
            elseif distance_cc(length_cc)<sqrt(3)*a+eps*1e2
                cc_neighbor_list(exh_b)=2;
            elseif distance_cc(length_cc) < 2*a+eps*1e2
                cc_neighbor_list(exh_b)=3;
            end
        end
        
    end
    n_dlt_ssgraph=find(dlt_ssgraph==1);
arrangement_ssgraph=arrangement_ssgraph(n_dlt_ssgraph,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%solve the least squares function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length_ssgraphcrd=size(arrangement_ssgraph,1);
trans_x=zeros(1,length_ssgraphcrd);
trans_y=zeros(1,length_ssgraphcrd);

trans_theta=zeros(1,length_ssgraphcrd);

for exh=1:length_ssgraphcrd
    fprintf('exh=%d\n',exh);
    ssgraph(1,:)=scgraph(1,arrangement_ssgraph(exh,:));
    ssgraph(2,:)=scgraph(2,arrangement_ssgraph(exh,:));
    
    trans_theta_1=[cc(1,2)-cc(1,1),cc(2,2)-cc(2,1)];
    trans_theta_0=[ssgraph(1,2)-ssgraph(1,1),ssgraph(2,2)-ssgraph(2,1)];
    trans_theta(exh)=get_angle(trans_theta_0,trans_theta_1);
    Trans_T=[cos(-trans_theta(exh)),-sin(-trans_theta(exh)),0;sin(-trans_theta(exh)),cos(-trans_theta(exh)),0;0,0,1];
    cc_temp1=Trans_T*cc;
    trans_x(exh)=mean(ssgraph(1,:))-mean(cc(1,:));
    trans_y(exh)=mean(ssgraph(2,:))-mean(cc(2,:));

    cc_temp=[cc_temp1(1,:)+trans_x(exh);cc_temp1(2,:)+trans_y(exh);cc_temp1(3,:)];
    if switch_structure_plot==1
        figure (500+exh)
        plot(ssgraph(1,:),ssgraph(2,:),'-*')
        hold on
        plot(cc(1,:),cc(2,:),'-s')
        hold on
        plot(cc_temp(1,:),cc_temp(2,:),'-o')
        hold on
        plot(cc_temp1(1,:),cc_temp1(2,:),'-^')
        hold off
        axis equal
        pause
    end
   
    [subf,X,Y,alpha1,beta1,gamma1,h]=JoiningCNT2(ssgraph,cc_temp,a,n,m,z);
    Xarray(exh)=X;
    Yarray(exh)=Y;
    alphaarray(exh)=alpha1;
    betaarray(exh)=beta1;
    gammaarray(exh)=gamma1;
    harray(exh)=h;
    subfarray(exh)=subf(1);
end

if length(subfarray)>=20
    subfarray_idx=find(subfarray>=0 & subfarray<mean(subfarray) & harray>0);
else
    subfarray_idx=find(subfarray>=0 & harray>0);
end
delflag_subfarray_idx=zeros(1,length(subfarray_idx));
num_overlap=zeros(1,length(subfarray_idx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Rank the results by bonding conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for exh=1:length(subfarray_idx)
    ssgraph(1,:)=scgraph(1,arrangement_ssgraph(subfarray_idx(exh),:));
    ssgraph(2,:)=scgraph(2,arrangement_ssgraph(subfarray_idx(exh),:));
    X=Xarray(subfarray_idx(exh));
    Y=Yarray(subfarray_idx(exh));
    h=harray(subfarray_idx(exh));
    alpha1=alphaarray(subfarray_idx(exh));
    beta1=betaarray(subfarray_idx(exh));
    gamma1=gammaarray(subfarray_idx(exh));
    ssgraph=[ssgraph(1,:)+X;ssgraph(2,:)+Y];
    
    scgraph_new_temp=[scgraph(1,:)+X;scgraph(2,:)+Y;zeros(1,length(scgraph))];
    scgraphnewx=scgraph_new_temp(1,:);
    scgraphnewy=scgraph_new_temp(2,:);
    [in_cnt_temp,on_cnt_temp] = inpolygon(scgraphnewx,scgraphnewy,ssgraph(1,:),ssgraph(2,:));
    del_cnt_temp=in_cnt_temp~=on_cnt_temp;
    scgraph_new_temp(:,del_cnt_temp)=[];
    
    Trans_T=[cos(-trans_theta(subfarray_idx(exh))),-sin(-trans_theta(subfarray_idx(exh))),0;sin(-trans_theta(subfarray_idx(exh))),cos(-trans_theta(subfarray_idx(exh))),0;0,0,1];
    cc_temp1=Trans_T*cc;
    cc_temp=[cc_temp1(1,:)+trans_x(subfarray_idx(exh));cc_temp1(2,:)+trans_y(subfarray_idx(exh));cc_temp1(3,:)];
    

    Trans_3D=[cos(alpha1)*cos(gamma1)-cos(beta1)*sin(alpha1)*sin(gamma1),-cos(beta1)*cos(gamma1)*sin(alpha1)-cos(alpha1)*sin(gamma1),sin(alpha1)*sin(beta1);
    cos(gamma1)*sin(alpha1)+cos(alpha1)*cos(beta1)*sin(gamma1),cos(alpha1)*cos(beta1)*cos(gamma1)-sin(alpha1)*sin(gamma1),-cos(alpha1)*sin(beta1);
    sin(beta1)*sin(gamma1), cos(gamma1)*sin(beta1),cos(beta1)
    ];
    
    cc_temp1=Trans_3D*cc_temp+[0;0;h];
    cc_new=cc_temp1(1:2,:);
    
    
    D_cnt_ssgrph=zeros(1,sizecc);
    
    [in_cnt_pjt,on_cnt_pjt] = inpolygon(cc_new(1,:),cc_new(2,:),ssgraph(1,:),ssgraph(2,:));
    %% number of carbon nanotube bottoms atoms inside the polygon formed by selected graphene atoms
    if sum(in_cnt_pjt)>=floor(sizecc*tolerance)
        delflag_subfarray_idx(exh)=0;
    elseif sum(in_cnt_pjt)<floor(sizecc*tolerance)
        if sum(in_cnt_pjt)<=0.5*floor(sizecc*tolerance)
            delflag_subfarray_idx(exh)=2;
        else
            delflag_subfarray_idx(exh)=1;
        end
    end
    if switch_reslt_plot==1
        figure(exh+20)%% plot all possible junctions
        axis equal
        plot(scgraph_new_temp(1,:),scgraph_new_temp(2,:),'k.','MarkerSize',30);
        hold on
        plot(scgraphnewx(del_cnt_temp),scgraphnewy(del_cnt_temp),'g.','MarkerSize',30);
        hold on
        plot(ssgraph(1,:),ssgraph(2,:),'b.','MarkerSize',30);
        hold on
        plot(cc_new(1,:),cc_new(2,:),'r.','MarkerSize',30)
        hold off
        pause
    end
    
    if delflag_subfarray_idx(exh)==0
        continue
    end
    for exh1=1:sizecc
        D_cnt_ssgrph(1,exh1)=norm(cc_new(:,exh1)-ssgraph(:,exh1));
    end
    
    D_cnt_grph=zeros(sizecc,length(scgraph_new_temp));
    for exh2=1:sizecc
        D_cnt_grph(exh2,:)=sqrt((scgraph_new_temp(1,:)-cc_new(1,exh2)).^2 + (scgraph_new_temp(2,:)-cc_new(2,exh2)).^2);
    end
    dlt_D_cnt_grph=find(in_cnt_temp==1);
        %% make sure no sp3 bonds by calculating the distance between carbon nanotube atoms and other graphene atoms
    for exh3=1:sizecc
        if sum(D_cnt_grph(exh3,:)<D_cnt_ssgrph(1,exh3)-1e2*eps)>1
            delflag_subfarray_idx(exh)=2;
        end
    end
    
    num_overlap(exh)=length(find( sqrt((ssgraph(1,:)-cc_new(1,:)).^2+(ssgraph(2,:)-cc_new(2,:)).^2)<0.1));
    %% calculate the average distance between each group of graphene and carbon nanotube atoms
    if (num_overlap(exh)>=sizecc_partial)&&(delflag_subfarray_idx(exh)<2)
        delflag_subfarray_idx(exh)=0;
    end
end

figure(6) %% plot the residual values
plot(subfarray(subfarray_idx),'o-')
for i = 1:length(subfarray_idx)
    text(i+0.05, subfarray(subfarray_idx(i)), num2str(subfarray(subfarray_idx(i))), 'tag', 'text');
end

if find(delflag_subfarray_idx==0)
    column0=find(delflag_subfarray_idx==0);
    minsubf=min(subfarray(subfarray_idx(column0)));
elseif find(delflag_subfarray_idx==1)
    column0=find(delflag_subfarray_idx==1);
    minsubf=min(subfarray(subfarray_idx(column0)));
    disp('Warning! The structure may be not stable enough!')
else
    column0=find(delflag_subfarray_idx==2);
    minsubf=min(subfarray(subfarray_idx(column0)));
    disp('Warning! The structure is very unstable!')
end

%%Determine the CNT/graphene junction
[row,column]=find(subfarray==minsubf);
X=Xarray(column(1));
Y=Yarray(column(1));
h=harray(column(1));
alpha1=alphaarray(column(1));
beta1=betaarray(column(1));
gamma1=gammaarray(column(1));

fprintf('X=%d\nY=%d\nh=%d\nalpha=%d\nbeta=%d\ngamma=%d\nf=%f\n',X,Y,h,alpha1,beta1,gamma1,minsubf);

figure (11)%%plot the selected graphene atoms 
Trans_T=[cos(-trans_theta(column(1))),-sin(-trans_theta(column(1))),0;sin(-trans_theta(column(1))),cos(-trans_theta(column(1))),0;0,0,1];
c_temp1=Trans_T*c;
Trans_3D=[cos(alpha1)*cos(gamma1)-cos(beta1)*sin(alpha1)*sin(gamma1),-cos(beta1)*cos(gamma1)*sin(alpha1)-cos(alpha1)*sin(gamma1),sin(alpha1)*sin(beta1);
    cos(gamma1)*sin(alpha1)+cos(alpha1)*cos(beta1)*sin(gamma1),cos(alpha1)*cos(beta1)*cos(gamma1)-sin(alpha1)*sin(gamma1),-cos(alpha1)*sin(beta1);
    sin(beta1)*sin(gamma1), cos(gamma1)*sin(beta1),cos(beta1)
    ];
    

c_final=c_temp1+[trans_x((column(1)));trans_y((column(1))); 0];
ctrans=Trans_3D*c_final;
cnew=[ctrans(1,:);ctrans(2,:);ctrans(3,:)+h]; %%CNT coordinate
cc_temp1=Trans_T*cc;
cc_temp=[cc_temp1(1,:)+trans_x(column(1));cc_temp1(2,:)+trans_y(column(1));cc_temp1(3,:)];
cctrans1=cc_temp(1:2,:);
cctrans2=cctrans1(1:2,:);

ssgraph(1,:)=scgraph(1,arrangement_ssgraph(column(1),:));
ssgraph(2,:)=scgraph(2,arrangement_ssgraph(column(1),:));
scgraphnew=[scgraph(1,:)+X;scgraph(2,:)+Y;zeros(1,length(scgraph))];
scgraphnewx=scgraphnew(1,:);
scgraphnewy=scgraphnew(2,:);
scgraphnewz=scgraphnew(3,:);

%% Lift graphene atoms if necessary
ssgraph=[ssgraph(1,:)+X;ssgraph(2,:)+Y];
[in_cnt,on_cnt] = inpolygon(scgraphnewx,scgraphnewy,ssgraph(1,:),ssgraph(2,:));
del_cnt=in_cnt~=on_cnt;
lift_atoms_cnt=zeros(1,sizecc);
neighbor_list_copy(del_cnt,:)=0;
neighbor_list_copy(:,del_cnt)=0;

for i_lift_atoms_cnt=1:sizecc
    lift_atoms_cnt(i_lift_atoms_cnt)=find(c(1,:)==cc(1,i_lift_atoms_cnt)&c(2,:)==cc(2,i_lift_atoms_cnt)&c(3,:)==cc(3,i_lift_atoms_cnt));
end
ccnew=Trans_3D*cc_temp+[0;0;h];
lift_atoms_dists=zeros(1,sizecc);
lift_atoms_dists_p=zeros(1,sizecc);
delta_lift_atoms_dists=zeros(1,sizecc);
slct_G=arrangement_ssgraph(column(1),:);

if (m~=n) && (m>=20)
    for i_lift_atoms_dists=1:sizecc
        lift_atoms_dists(i_lift_atoms_dists)=norm([ssgraph(:,i_lift_atoms_dists);0]-ccnew(:,i_lift_atoms_dists));  %%L
        lift_atoms_dists_p(i_lift_atoms_dists)=norm(ssgraph(:,i_lift_atoms_dists)-ccnew(1:2,i_lift_atoms_dists));  %%t
        
        if lift_atoms_dists(i_lift_atoms_dists) < D_cutoff
            delta_lift_atoms_dists(:,i_lift_atoms_dists)=abs(sqrt(lift_atoms_dists(i_lift_atoms_dists)^2- lift_atoms_dists_p(i_lift_atoms_dists)^2))- ...
                abs(sqrt((a)^2- lift_atoms_dists_p(i_lift_atoms_dists)^2));
        else
            nn_lift_atoms_dists=find(neighbor_list_copy(slct_G(i_lift_atoms_dists),:)==1);
            delta_lift_atoms_dists(:,i_lift_atoms_dists)=ccnew(3,i_lift_atoms_dists)/3;
        end
    end
    
    ssgraph=[ssgraph;delta_lift_atoms_dists];
    scgraphnew(3,arrangement_ssgraph(column(1),:))=delta_lift_atoms_dists;
    lift_atoms_dists_after=zeros(1,sizecc);
    for i_lift_atoms_dists=1:sizecc
        lift_atoms_dists_after(i_lift_atoms_dists)=norm(ssgraph(:,i_lift_atoms_dists)-ccnew(1:3,i_lift_atoms_dists));
    end

end

plot(ssgraph(1,:),ssgraph(2,:),'ko')
hold off




figure(10) %%plot the deleted atoms on graphene
plot(scgraphnewx(in_cnt),scgraphnewy(in_cnt),'go')%%?????????
hold on
plot(scgraphnewx(del_cnt),scgraphnewy(del_cnt),'rx')%%?????????
hold off
scgraphnew(:,del_cnt)=[];

figure(7) %%plot the junction
plot3(scgraphnew(1,:),scgraphnew(2,:),scgraphnew(3,:),'k.','MarkerSize',30)
axis equal
hold on
plot(ssgraph(1,:),ssgraph(2,:),'b.','MarkerSize',30)
axis equal
hold on
plot3(ccnew(1,:),ccnew(2,:),ccnew(3,:),'r.','MarkerSize',30)
hold off

%%output junction file for ovito
newjctfile=['Junction(',num2str(n),',',num2str(m),',',num2str(z),').dat'];
delete(newjctfile);

fjctid=fopen(newjctfile,'a');
leng1=length(scgraphnew);
leng2=length(ccnew);
leng3=leng1+leng2;
Out1=[num2str(leng3),' atoms\n'];
Out2=[num2str(min(scgraphnew(1,:))-1.5),' ',num2str(max(scgraphnew(1,:))+1.5),' xlo xhi\n'];
Out3=[num2str(min(scgraphnew(2,:))-1.5),' ',num2str(max(scgraphnew(2,:))+1.5),' ylo yhi\n'];
Out4=[num2str(min(scgraphnew(3,:))-1.5),' ',num2str(max(scgraphnew(3,:))+1.5),' zlo zhi\n'];
fprintf(fjctid,'# LAMMPS data file\n');
fprintf(fjctid,Out1);
fprintf(fjctid,'2 atom types\n');
fprintf(fjctid,Out2);
fprintf(fjctid,Out3);
fprintf(fjctid,Out4);
fprintf(fjctid,'\n');
fprintf(fjctid,'Atoms\n');
fprintf(fjctid,'\n');
atom1=1;
atom2=2;
for j1=1:1:leng1
    fprintf(fjctid,'%d %d %12.6f %12.6f %12.6f\n', j1, atom1, scgraphnew(1,j1), scgraphnew(2,j1), scgraphnew(3,j1));
end
for j2=1:1:leng2
    fprintf(fjctid,'%d %d %12.6f %12.6f %12.6f\n', leng1+j2, atom2, ccnew(1,j2), ccnew(2,j2), ccnew(3,j2));
end
fclose(fjctid);

figure(12)%%plot final CNT/graphene sturcture (red CNT, black graphene)
plot3(scgraphnew(1,:),scgraphnew(2,:),scgraphnew(3,:),'k.','MarkerSize',30)
hold on
plot3(cnew(1,:),cnew(2,:),cnew(3,:),'r.','MarkerSize',30)
hold off
figure (9)%%plot final CNT/graphene sturcture
cntgrph=[scgraphnew(1,:),cnew(1,:);scgraphnew(2,:),cnew(2,:);scgraphnew(3,:),cnew(3,:)];
plot3(cntgrph(1,:),cntgrph(2,:),cntgrph(3,:),'ko','MarkerSize',5)

xxxx=cntgrph(1,:);
yyyy=cntgrph(2,:);
zzzz=cntgrph(3,:);

%%output final CNT/graphene sturcture for ovito
newfile=['CNTgrph(',num2str(n),',',num2str(m),',',num2str(z),').dat'];
delete(newfile);

fid=fopen(newfile,'a');
leng1=length(scgraphnew);
leng2=length(cnew);
leng=length(xxxx);
Out1=[num2str(leng),' atoms\n'];
Out2=[num2str(min(xxxx)-1.5),' ',num2str(max(xxxx)+1.5),' xlo xhi\n'];
Out3=[num2str(min(yyyy)-1.5),' ',num2str(max(yyyy)+1.5),' ylo yhi\n'];
Out4=[num2str(min(zzzz)-1.5),' ',num2str(max(zzzz)+1.5),' zlo zhi\n'];
fprintf(fid,'# LAMMPS data file\n');
fprintf(fid,Out1);
fprintf(fid,'2 atom types\n');
fprintf(fid,Out2);
fprintf(fid,Out3);
fprintf(fid,Out4);
fprintf(fid,'\n');
fprintf(fid,'Atoms\n');
fprintf(fid,'\n');
atom1=1;
atom2=2;
for j1=1:1:leng1
    fprintf(fjctid,'%d %d %12.6f %12.6f %12.6f\n', j1, atom1, scgraphnew(1,j1), scgraphnew(2,j1), scgraphnew(3,j1));
end
for j2=1:1:leng2
    fprintf(fjctid,'%d %d %12.6f %12.6f %12.6f\n', leng1+j2, atom2, cnew(1,j2), cnew(2,j2), cnew(3,j2));
end

fclose(fid);

%%%%%%%Output XYZ file of final CNT/Graphene structure
newxyzfile=['CNTgrph',num2str(n),num2str(m),num2str(z),'.xyz'];
delete(newxyzfile);

fxyzid=fopen(newxyzfile,'a');
leng=length(xxxx);
Out1=num2str(leng);
Out2=['CNTgrph',num2str(n),num2str(m),num2str(z),'.xyz'];

fprintf(fxyzid,Out1);
fprintf(fid,'\n');
fprintf(fxyzid,Out2);
fprintf(fid,'\n');
for j=1:1:leng
    fprintf(fxyzid,'C %12.6f %12.6f %12.6f\n', xxxx(j), yyyy(j), zzzz(j));
end

fclose(fxyzid);
molviewer(newxyzfile)
%%%%
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%generate the rectangle of graphene for carbon nanotube.
function [CD,x1,x2,y1,y2]=Rec(n,m,zz)
a=1.42;
sign=0;
y(1,1) = 0;
x(1,1) = a*cos(pi/6);
CD=zeros(1,2);
for i = 2:4
    if mod(i,2)==1
        y(1,i) = y(1,i-1) + a;
    else
        y(1,i) = y(1,i-1) + a/2;
        x(1,i) = sign;
        if (i+1)<= 4
            x(1,i+1) = sign;
        end
        if sign == 0
            sign =a*cos(pi/6);
        else
            sign = 0;
        end
    end
end
for i = 2: 2
    for j = 1: 4
        x(i,j) = x(i-1,j) + 2*a*cos(pi/6);
        y(i,j) = y(1,j);
    end
end

x=x-a*cos(pi/6);
n0=[x(1,3),y(1,3)];
m0=[x(2,3),y(2,3)];
CD=n*n0+m*m0;


if m==n
    x1=zz;
    y1= CD(1, 2);
    x2=zz;
    y2=y(1,1)+(y1-CD(1,2));
else
    x1=(CD(1, 1)^2 + CD(1, 2)^2 - (CD(1, 2) - zz*(1/(CD(1, 1)^2 + CD(1, 2)^2))^(1/2)*CD(1, 1))*CD(1, 2))/CD(1, 1);
    y1= CD(1, 2) - zz*(1/(CD(1, 1)^2 + CD(1, 2)^2))^(1/2)*CD(1, 1);
    x2=x(1,1)+(x1-CD(1,1));
    y2=y(1,1)+(y1-CD(1,2));
end
end


%% rotate the retangle and make a side on the X axis

function [xnew,ynew]=tl(xxx,yyy,CD,oy)
a = [xxx';yyy'];
n=acos(dot(oy,CD)/(norm(CD)*norm(oy)));
T=[cos(n),-sin(n);sin(n),cos(n)];
A=T*a;
xnew=A(1,:);
ynew=A(2,:);
figure (3); plot(xnew,ynew,'x')%% plot graphene
end

function [output,line]=guess_next(current_index,current_num,line,output)
global n m N_neighbor slct_atom_index;
if (current_num==n+m)
    if check_slct(1)
        output(line,:)=slct_atom_index;
        line=line+1;
        return;
    else
        return;
    end
end


for i_neighbor=1:N_neighbor(current_index)
    
    new_index= get_neighbor_index(current_index,i_neighbor);
    
    
    if check_index(new_index,current_index, current_num)
        
        new_num=current_num+1;
        slct_atom_index(new_num)=new_index;
        
        [output,line]=guess_next(new_index,new_num,line,output);
    end
end

end



function new_index=get_neighbor_index(current_index,i)
global neighbor_list
index_list=find(neighbor_list(current_index,:)~=0);
new_index=index_list(i);
end


function n_check_index=check_index(new_index,current_index, current_num)
global scgraph slct_atom_index n m a any_two_distance ssgraphatom_0

x_new=scgraph(1,new_index);
y_new=scgraph(2,new_index);
x_current=scgraph(1,current_index);
y_current=scgraph(2,current_index);

n_check_index=1;
%%check postion of first guess
if current_num==1
    if (y_new>=y_current && x_new<=x_current)%
        return;
    else
        n_check_index=0;
        return;
    end
end

if current_num>=ceil((n+m)/2)
    if any_two_distance(new_index,slct_atom_index(1))>ceil((n+m)*3/4)*a
        n_check_index=0;
        return
    end
end

if current_num==n+m-1
    if any_two_distance(new_index,slct_atom_index(1))>2*a+1e2*eps
        n_check_index=0;
        return
    end
    if (ssgraphatom_0(2)-1e2*eps<=y_current && ssgraphatom_0(2)+1e2*eps>=y_current && ssgraphatom_0(1)<=x_current)
        n_check_index=0;
        return
    end
end

if sum(slct_atom_index==new_index)>0
    n_check_index=0;
    return
end


%second angle check
vec0=[1,0];
x_1=scgraph(1,slct_atom_index(1))-a;
y_1=scgraph(2,slct_atom_index(1));
vec1=[x_current-x_1,y_current-y_1];
vec2=[x_new-x_1,y_new-y_1];
angle1=get_angle(vec0,vec1);
angle2=get_angle(vec0,vec2);
if angle2<angle1
    n_check_index=0;
    return
end

%first angle check
old_index=slct_atom_index(current_num-1);
x_old=scgraph(1,old_index);
y_old=scgraph(2,old_index);
vec_0=[x_current-x_old,y_current-y_old];
vec_1=[x_new-x_current,y_new-y_current];
angle=get_angle(vec_0,vec_1);

if ~(angle>=0 && angle<=150/180*pi) %%take care of the angle
    n_check_index=0;
    return;
end



return;

end

function angle=get_angle(vec_0,vec_1)
over_pi=cross([vec_0,0],[vec_1,0]);
if over_pi(3)>=0
    angle=acos(dot(vec_0,vec_1)/norm(vec_0)/norm(vec_1));
else
    angle=2*pi-acos(dot(vec_0,vec_1)/norm(vec_0)/norm(vec_1));
end

if ~isreal(angle)||isnan(angle)
    angle=0;
end
end

function n_check_slct=check_slct(n)
global slct_atom_index
%%%For more function
n_check_slct=1;

end
