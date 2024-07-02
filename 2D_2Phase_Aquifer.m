clc
clear
%%%%%%%%%%%%%%%%%%%%%input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=[0	509	491	596	0	0	0	0	0	0	0	0
    439	509	491	596	526	0	0	0	439	877	0	0
    439	509	491	596	526	561	912	807	439	877	544	772
    439	509	491	596	526	561	912	807	439	877	544	772
    439	509	491	596	526	561	912	807	439	877	0	0
    0	0	0	596	526	561	912	807	439	877	0	0
    0	0	0	0	0	561	912	807	0	0	0	0
    0	0	0	0	0	0	912	807	439	0	0	0
    0	0	0	0	0	0	912	807	439	0	0	0];
dy=[0	474	474	474	0	0	0	0	0	0	0	0
    404	404	404	404	404	0	0	0	404	404	0	0
    386	386	386	386	386	386	386	386	386	386	386	386
    491	491	491	491	491	491	491	491	491	491	491	491
    404	404	404	404	404	404	404	404	404	404	0	0
    0	0	0	316	316	316	316	316	316	316	0	0
    0	0	0	0	0	316	316	316	0	0	0	0
    0	0	0	0	0	0	421	421	421	0	0	0
    0	0	0	0	0	0	526	526	526	0	0	0];
depth=[0	9342 9345 9347 0	0	 0    0	   0    0    0	  0
       9341	9327 9330 9338 9333	0	 0    0	   9311	9310 0	  0
       9336	9319 9316 9322 9325	9315 9299 9300 9299	9297 9297 9305
       9340	9326 9316 9308 9310	9323 9297 9296 9295	9293 9292 9295
       9342	9332 9323 9305 9292	9298 9296 9292 9291	9288 0	  0
       0	0	 0	  9315 9297	9295 9292 9289 9289	9287 0	  0
       0	0	 0	  0	   0	9294 9290 9286 0    0	 0	  0
       0	0  	 0	  0	   0	0	 9289 9281 9282	0	 0    0
       0	0	 0	  0	   0	0	 9290 9280 9278	0    0    0];
dz=[0	10	12	5	0	0	0	0	0	0	0	0
    8	35	30	15	6	0	0	0	4	5	0	0
    14	44	36	30	22	16	12	14	15	11	6	3
    20	34	35	40	34	32	29	25	22	18	10	3
    5	12	12	40	44	42	32	20	16	10	0	0
    0	0	0	10	19	27	24	10	6	3	0	0
    0	0	0	0	0	4	10	6	0	0	0	0
    0	0	0	0	0	0	8	7	3	0	0	0
    0	0	0	0	0	0	4	5	2	0	0	0];
kx=[0	275	270	252	0	0	0	0	0	0	0	0
    267	274	280	265	253	0	0	0	259	270	0	0
    265	280	289	278	271	271	270	269	270	279	283	275
    258	271	295	297	282	280	281	276	290	293	279	270
    253	259	275	285	290	280	289	277	290	280	0	0
    0	0	0	272	276	273	288	281	274	268	0	0
    0	0	0	0	0	265	280	290	0	0	0	0
    0	0	0	0	0	0	270	280	270	0	0	0
    0	0	0	0	0	0	260	268	260	0	0	0];
ky=[0	  220	216   201.6 0	  0	    0     0	    0	  0	    0	  0
    213.6 219.2	224	  212	202.4 0	    0	  0	    207.2 216	0	  0
    212	  224	231.2 222.4	216.8 216.8	216	  215.2	216	  223.2	226.4 220
    206.4 216.8	236	  237	225.6 224	224	  220.8	232	  234.4	223.2 216
    202.4 207.2	220	  228	232	  224	231.2 221.6	232	  224	0	  0
    0	  0	    0	  217.6	220.8 218.4	230.4 224.8	219.2 214.4	0	  0
    0	  0	    0	  0	    0	  212	224	  232	0	  0	    0	  0
    0	  0	    0	  0	    0	  0	    216	  224	216	  0	    0	  0
    0	  0	    0	  0	    0	  0	    208	  214.4	208	  0	    0	  0];
por=[0	    0.192	0.197	0.202	0	    0	    0	    0    	0	    0	    0	    0
     0.19	0.195	0.2	    0.204	0.207	0	    0	    0	    0.215	0.205	0	    0
     0.19	0.196	0.205	0.207	0.21	0.216	0.22	0.223	0.215	0.21	0.203	0.2
     0.185	0.195	0.205	0.213	0.216	0.221	0.225	0.226	0.22	0.215	0.207	0.2
     0.183	0.195	0.205	0.212	0.218	0.225	0.232	0.232	0.225	0.219	0	    0
     0	    0	    0	    0.21	0.219	0.226	0.235	0.23	0.22	0.216	0	    0
     0	    0	    0     	0	    0	    0.225	0.235	0.23	0	    0	    0	    0
     0	    0	    0	    0	    0	    0	    0.232	0.226	0.217	0	    0	    0
     0	    0     	0	    0    	0	    0	    0.229	0.22	0.217	0	    0	    0];
%%%%initial oil pressure%%%%if aquifer -1 and if no flow B.C 0%%%%%%%%%%%%%
poi=[0	0	 0	  0	   0	0	 0	  0    0	0	 0	  0	   0	0
     0	0	 7000 7000 7000	0	 0	  0	   0	0	 0	  0	   0	0
     0	7000 7000 7000 7000	7000 0	  0	   0	7000 7000 0	   0	0
     0	7000 7000 7000 7000	7000 7000 7000 7000	7000 7000 7000 7000	0
     0	7000 7000 7000 7000	7000 7000 7000 7000	7000 7000 7000 7000	-1
     -1	7000 7000 7000 7000	7000 7000 7000 7000	7000 7000 -1  -1	1
     1	-1	 -1	  -1   7000	7000 7000 7000 7000	7000 7000 -1   1	1
     1	1	 1	  1	   -1	-1	 7000 7000 7000	-1	 -1	  1	   1	1
     1	1	 1	  1	   1	1	 -1	  7000 7000 7000 -1	  1	   1	1
     1	1	 1	  1	   1	1	 -1	  7000 7000 7000 -1	  1	   1	1
     1	1	 1	  1	   1	1	 1	  -1    -1	-1	 1	  1	   1	1]
 soi=[0.0	0.5	0.5	0.5	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
      0.5	0.5	0.5	0.5	0.5	0.0	0.0	0.0	0.5	0.5	0.0	0.0
      0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5
      0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5
      0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.0	0.0
      0.0	0.0	0.0	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.0	0.0
      0.0	0.0	0.0	0.0	0.0	0.5	0.5	0.5	0.0	0.0	0.0	0.0
      0.0	0.0	0.0	0.0	0.0	0.0	0.5	0.5	0.5	0.0	0.0	0.0
      0.0	0.0	0.0	0.0	0.0	0.0	0.5	0.5	0.5	0.0	0.0	0.0];
pc=[.2 8
    .25 4.3
    .3 3
    .4 1.78
    .5 1.21
    .6 .79
    .7 .43
    .8 .1
    .9 0];
oilvis=[5000 .92
    5500 .9243
    6000 .9372
    6500 .9494
    7000 .965
    7500 .9812
    8000 1.0019];
watvis=[5000 .52
        6000 .52];
sat=[.18 0 1
    .21 .00000 .92692
    .24 .00002 .85441
    .27 .00014 .79288
    .30 .00045 .71312
    .33 .00111 .64526
    .36 .00232 .57980
    .39 .00430 .51709
    .42 .00733 .45744
    .45 .01175 .40110
    .48 .01791 .34831
    .51 .02623 .29924
    .54 .03714 .25403
    .57 .05116 .21278
    .60 .06882 .17552
    .63 .09069 .14228
    .66 .11741 .11301
    .69 .14963 .08763
    .72 .18807 .06603
    .75 .23347 .04803
    .78 .28664 .03344
    .81 .34842 .02199
    .84 .41968 .01340
    .87 .50135 .00733
    .90 .59439 .00340];
 %%%%%%%%%%main dimention of the grid system i:downward j:forwrd%%%%%%%%%%
dim=[9 12];
dens=[45 62];
cr=[5.9e-6];
%well=[icoord,jcoord,oil=1 wat=0,ratecte=1 pwfcte=0,rw,pwf,rate prod<0 inj>0;..]
well=[2 3 1 0 .25 5300 -1;3 9 1 0 .25 5300 -1;4 6 0 1 .25 1 2500];
%%%%%%%%%%%%%%%%%%%delta t and final time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=[5 5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%grid center depth%%%%%%%%%%%%%%%%%%%
sizepc=size(pc);
voi=0;
for i=1:dim(1)
    for j=1:dim(2)
        z(i,j)=depth(i,j)+0.5*dz(i,j);
        ii=i+1;
        jj=j+1;
        poil(i,j)=poi(ii,jj);
        sato(i,j)=soi(i,j);
        satw(i,j)=1-sato(i,j);
        pcow(i,j)=0;
        
        if satw(i,j)>pc(sizepc(1),1)
            pcow(i,j)=0;
        else
        pcow(i,j)=interp1(pc(:,1),pc(:,2),satw(i,j),'linear'); 
        end
        
        pwat(i,j)=poil(i,j)-pcow(i,j);
        
        if dx(i,j)==0
    poil(i,j)=0;pwat(i,j)=0;
        end
        
        boil(i,j)=(1/(1+5e-6*(poil(i,j)-14.7)));
        
        voi=(((dx(i,j).*dy(i,j).*dz(i,j))./5.615).*por(i,j).*(sato(i,j)./boil(i,j)))+voi;
    end
end
rec=0;
step=0;
sizesat=size(sat);
qocum=0;
qwcum=0;
qaqucum=0;
qaqu=zeros(dim(1),dim(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%main loop%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              for t=time(1):time(1):time(2)
error=11;
poiln=poil-100;
pcown=pcow;
step=step+1
                   while abs(error)>9

err(i,j)=0;
erro=0;
error;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pysical properties matrix
for i=1:dim(1)
  for j=1:dim(2)
%%%%%%%%%%%%%%%oil viscosity %%%%%water viscosity%%%%%%%%%%%%%%  
zmo=polyfit(oilvis(:,1),oilvis(:,2),1);
moil(i,j)=zmo(1,1).*poil(i,j)+zmo(1,2);
zmw=polyfit(watvis(:,1),watvis(:,2),1);
mwat(i,j)=zmw(1,1).*pwat(i,j)+zmw(1,2);
%%%%%%%%%%%%%%%oil FVF %%%%%%%water FVF %%%%%%%%%%%%%%%%%%%%%%
boil(i,j)=(1/(1+5e-6*(poil(i,j)-14.7)));
boiln(i,j)=(1/(1+5e-6*(poiln(i,j)-14.7)));
pwatn(i,j)=poiln(i,j)-pcown(i,j);
bwat(i,j)=(1/(1+1e-6*(pwat(i,j)-14.7)));
bwatn(i,j)=(1/(1+1e-6*(pwatn(i,j)-14.7)));
porn(i,j)=por(i,j).*(1+cr(1).*((poiln(i,j)+pwatn(i,j))./2-((poil(i,j)+pwat(i,j))./2)));
%%%%%%%%%%%%%%%%%kr oil%%%%%%%%%kr water%%%%%%%%%%%%%%%%%%%%%
if satw(i,j)<=sat(1,1)
    kroil(i,j)=1;
    krwat(i,j)=0;
elseif satw(i,j)>=sat(sizesat(1),1)
    kroil(i,j)=sat(sizesat(1),3);
    krwat(i,j)=sat(sizesat(1),2);
else
   kroil(i,j)=interp1(sat(:,1),sat(:,3),satw(i,j),'linear');
   krwat(i,j)=interp1(sat(:,1),sat(:,2),satw(i,j),'linear');
end  
if dx(i,j)==0
    kroil(i,j)=0;krwat(i,j)=0;
    bwat(i,j)=0;boil(i,j)=0;
    mwat(i,j)=0;moil(i,j)=0;
    poil(i,j)=0;pwat(i,j)=0;
    bwatn(i,j)=0;boiln(i,j)=0;
    satw(i,j)=0;sato(i,j)=0;
end
%%%%%%%%%%%%%%%landao landaw landasumogw%%%%(oil water mobility)%%%%%%%
if or(boil(i,j)==0,moil(i,j)==0)
    landao(i,j)=0;
else
landao(i,j)=(kroil(i,j)./(boil(i,j).*moil(i,j)));
end
if or(bwat(i,j)==0,mwat(i,j)==0)
    landaw(i,j)=0;
else
landaw(i,j)=(krwat(i,j)./(bwat(i,j).*mwat(i,j)));
end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% well%%%%%%%%%%%%%%%%%%%%
%well=[icoord,jcoord,oil=1 wat=0,ratecte=1 pwfcte=0,rw,pwf,rate prod<0 inj>0;...]
qo=zeros(dim(1),dim(2));
qw=zeros(dim(1),dim(2));
sizewell=size(well);
for kk=1:sizewell(1)
%%%%%%%%%%%%%%PWF CTE
    if well(kk,4)==0
        if or(kx(well(kk,1),well(kk,2))==0,ky(well(kk,1),well(kk,2))==0)
           pid(well(kk,1),well(kk,2))=0;
        else
        requ(well(kk,1),well(kk,2))=(.28*(((ky(well(kk,1),well(kk,2))./kx(well(kk,1),well(kk,2))).^.5).*(dx(well(kk,1),well(kk,2))).^2+...
                                              ((kx(well(kk,1),well(kk,2))./ky(well(kk,1),well(kk,2))).^.5).*(dy(well(kk,1),well(kk,2))).^2)^.5)./...
                                             (((ky(well(kk,1),well(kk,2))./kx(well(kk,1),well(kk,2))).^.25)+...
                                              ((kx(well(kk,1),well(kk,2))./ky(well(kk,1),well(kk,2))).^.25));
         kabs(well(kk,1),well(kk,2))=((ky(well(kk,1),well(kk,2)).*kx(well(kk,1),well(kk,2))).^.5);                           
         pid(well(kk,1),well(kk,2))=(0.0078.*kabs(well(kk,1),well(kk,2)).*dz(well(kk,1),well(kk,2)))./(reallog(requ(well(kk,1),well(kk,2))./well(kk,5)));     
        end
      if well(kk,3)>-1                 % oil & water     
         if well(kk,7)<0       % oil & water production  
            qo(well(kk,1),well(kk,2))=-landao(well(kk,1),well(kk,2)).*(poil(well(kk,1),well(kk,2))-well(kk,6)).*pid(well(kk,1),well(kk,2));
            qw(well(kk,1),well(kk,2))=-landaw(well(kk,1),well(kk,2)).*(pwat(well(kk,1),well(kk,2))-well(kk,6)).*pid(well(kk,1),well(kk,2));
         elseif well(kk,7)>0       %    injection
             if well(kk,3)==1     % oil injection 
                 qo(well(kk,1),well(kk,2))=pid(well(kk,1),well(kk,2)).*(poil(well(kk,1),well(kk,2))-well(kk,6)).*...
                     (landao(well(kk,1),well(kk,2))+landaw(well(kk,1),well(kk,2)).*(bwat(well(kk,1),well(kk,2))./boil(well(kk,1),well(kk,2))));
             elseif well(kk,3)==0   %water injection  
                 qw(well(kk,1),well(kk,2))=pid(well(kk,1),well(kk,2)).*(pwat(well(kk,1),well(kk,2))-well(kk,6)).*...
                     (landaw(well(kk,1),well(kk,2))+landao(well(kk,1),well(kk,2)).*(boil(well(kk,1),well(kk,2))./bwat(well(kk,1),well(kk,2))));
             end
         end
      end
%%%%%%RATE CTE%%%%%%        
    elseif well(kk,4)==1
     if well(kk,3)==1                 % oil
        if well(kk,7)<0       % oil production  
          if landao(well(kk,2),well(kk,1))==0 
          continue
          end 
           qo(well(kk,1),well(kk,2))=well(kk,7);
           qw(well(kk,1),well(kk,2))=qo(well(kk,1),well(kk,2)).*(landaw(well(kk,1),well(kk,2))./landao(well(kk,1),well(kk,2)));                             
        elseif well(kk,7)>0    %   oil injection
          qo(well(kk,1),well(kk,2))=well(kk,7);           
        end                    
     elseif well(kk,3)==0           %   water
        if well(kk,7)<0      %   water production
           if landaw(well(kk,1),well(kk,2))==0 
           continue
           end
           qw(well(kk,1),well(kk,2))=well(kk,7);
           qo(well(kk,1),well(kk,2))=qo(well(kk,1),well(kk,2)).*(landao(well(kk,1),well(kk,2))./landaw(well(kk,1),well(kk,2)));                      
        elseif well(kk,7)>0   %  water injection
          qw(well(kk,1),well(kk,2))=well(kk,7); 
          requ(well(kk,1),well(kk,2))=(.28*(((ky(well(kk,1),well(kk,2))./kx(well(kk,1),well(kk,2))).^.5).*(dx(well(kk,1),well(kk,2))).^2+...
                                              ((kx(well(kk,1),well(kk,2))./ky(well(kk,1),well(kk,2))).^.5).*(dy(well(kk,1),well(kk,2))).^2)^.5)./...
                                             (((ky(well(kk,1),well(kk,2))./kx(well(kk,1),well(kk,2))).^.25)+...
                                              ((kx(well(kk,1),well(kk,2))./ky(well(kk,1),well(kk,2))).^.25));
         kabs(well(kk,1),well(kk,2))=((ky(well(kk,1),well(kk,2)).*kx(well(kk,1),well(kk,2))).^.5);                           
         pid(well(kk,1),well(kk,2))=(0.0078.*kabs(well(kk,1),well(kk,2)).*dz(well(kk,1),well(kk,2)))./(reallog(requ(well(kk,1),well(kk,2))./well(kk,5)));     
         pwf=pwat(well(kk,1),well(kk,2))-(qw(well(kk,1),well(kk,2))./(-landaw(well(kk,1),well(kk,2)).*pid(well(kk,1),well(kk,2))));
        end   
     end  
%%%%%%%%ERROR%%%%%   
    elseif well(kk,4)==-1
        disp('"error in well input"');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%transsmisibility calculation%%%%%%%%%%%%%
for i=1:dim(1)
  for j=1:dim(2)
       if dx(i,j)==0
        continue
       end
    if poi(i,j+1)<=0
        ip=i;
    else
        ip=i-1;
    end
    if poi(i+2,j+1)<=0
        in=i;
    else
        in=i+1;
    end
    if poi(i+1,j)<=0
        jp=j;
    else
        jp=j-1;
    end
    if poi(i+1,j+2)<=0
        jn=j;
    else
        jn=j+1;
    end
       toxp(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(ip,j).*dy(ip,j).*kx(ip,j))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(ip,j)+dz(ip,j).*dy(ip,j).*kx(ip,j).*dx(i,j))).*...                                         
                     (((kroil(i,j)+kroil(ip,j))./2)/(((boil(i,j)+boil(ip,j))./2).*(moil(i,j)+moil(ip,j))./2));
       twxp(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(ip,j).*dy(ip,j).*kx(ip,j))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(ip,j)+dz(ip,j).*dy(ip,j).*kx(ip,j).*dx(i,j))).*...                                         
                     (((krwat(i,j)+krwat(ip,j))./2)/(((bwat(i,j)+bwat(ip,j))./2).*(mwat(i,j)+mwat(ip,j))./2));
       toxn(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(in,j).*dy(in,j).*kx(in,j))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(in,j)+dz(in,j).*dy(in,j).*kx(in,j).*dx(i,j))).*...                                         
                     (((kroil(i,j)+kroil(in,j))./2)/(((boil(i,j)+boil(in,j))./2).*(moil(i,j)+moil(in,j))./2));
       twxn(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(in,j).*dy(in,j).*kx(in,j))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(in,j)+dz(in,j).*dy(in,j).*kx(in,j).*dx(i,j))).*...                                         
                     (((krwat(i,j)+krwat(in,j))./2)/(((bwat(i,j)+bwat(in,j))./2).*(mwat(i,j)+mwat(in,j))./2));
       toyp(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(i,jp).*dy(i,jp).*kx(i,jp))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(i,jp)+dz(i,jp).*dy(i,jp).*kx(i,jp).*dx(i,j))).*...                                         
                     (((kroil(i,j)+kroil(i,jp))./2)/(((boil(i,j)+boil(i,jp))./2).*(moil(i,j)+moil(i,jp))./2));
       twyp(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(i,jp).*dy(i,jp).*kx(i,jp))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(i,jp)+dz(i,jp).*dy(i,jp).*kx(i,jp).*dx(i,j))).*...                                         
                     (((krwat(i,j)+krwat(i,jp))./2)/(((bwat(i,j)+bwat(i,jp))./2).*(mwat(i,j)+mwat(i,jp))./2));
       toyn(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(i,jn).*dy(i,jn).*kx(i,jn))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(i,jn)+dz(i,jn).*dy(i,jn).*kx(i,jn).*dx(i,j))).*...                                         
                     (((kroil(i,j)+kroil(i,jn))./2)/(((boil(i,j)+boil(i,jn))./2).*(moil(i,j)+moil(i,jn))./2));
       twyn(i,j)=1.127.*2e-3.*((dz(i,j).*dy(i,j).*kx(i,j).*dz(i,jn).*dy(i,jn).*kx(i,jn))./...
                     (dz(i,j).*dy(i,j).*kx(i,j).*dx(i,jn)+dz(i,jn).*dy(i,jn).*kx(i,jn).*dx(i,j))).*...                                         
                     (((krwat(i,j)+krwat(i,jn))./2)/(((bwat(i,j)+bwat(i,jn))./2).*(mwat(i,j)+mwat(i,jn))./2));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%cop cow cwp cww calculation%%%%%%%%%%%%%
 for i=1:dim(1)
  for j=1:dim(2) 
    if dx(i,j)==0
        continue
    end
dpor(i,j)=(porn(i,j)-por(i,j))./(poiln(i,j)-poil(i,j));
dboil(i,j)=((1./boiln(i,j))-(1./boil(i,j)))./(poiln(i,j)-poil(i,j));
dbwat(i,j)=((1./bwatn(i,j))-(1./bwat(i,j)))./(poiln(i,j)-poil(i,j));
cop(i,j)=((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*(1-satw(i,j)).*...
    ((dpor(i,j)./boil(i,j))+porn(i,j).*dboil(i,j));
cow(i,j)=-((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*(porn(i,j)./boiln(i,j));
cwp(i,j)=((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*satw(i,j).*...
    ((dpor(i,j)./bwat(i,j))+porn(i,j).*dbwat(i,j));
cww(i,j)=((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*(porn(i,j)./bwatn(i,j));
  end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%poiln calculation%%%%%%%%%%%%%%%%%%%%
counter=0;
f=1;
u=1;
w=zeros(dim(1),dim(2));
zip=zeros(dim(1),dim(2));
zin=zeros(dim(1),dim(2));
zjn=zeros(dim(1),dim(2));
zjp=zeros(dim(1),dim(2));
for j=1:dim(2)
   for i=1:dim(1)
       w(i,j)=f-1;
    if dx(i,j)==0
        w(i,j)=f;
        f=f+1;
        continue
    end 
   end
end
for j=1:dim(2)
   for i=1:dim(1)
       i;j;
    if dx(i,j)==0
        u=u+1;
        continue
    end 
counter=counter+1;
position=i+(j-1)*dim(1)-u+1;
iip=(i-1)+(j-1)*dim(1);
jjp=i+(j-2)*dim(1);
jjn=i+(j)*dim(1);
iin=i+1+(j-1)*dim(1);
                 %%%%%%%%%%%%%%%%matrise zarayeb%%%%%%%%%%
           if i~=dim(1)
               if dx(i+1,j)~=0
       zin(i,j)=boiln(i,j).*toxn(i,j)+bwatn(i,j).*twxn(i,j);  
           a(counter,iin-w(i+1,j))=zin(i,j);
               end
           end   
           if i~=1
               if dx(i-1,j)~=0
       zip(i,j)=boiln(i,j).*toxp(i,j)+bwatn(i,j).*twxp(i,j);
           a(counter,iip-w(i-1,j))=zip(i,j);
               end
           end       
           if j~=1
               if dx(i,j-1)~=0
       zjp(i,j)=boiln(i,j).*toyp(i,j)+bwatn(i,j).*twyp(i,j); 
           a(counter,jjp-w(i,j-1))=zjp(i,j);
               end
           end     
           if j~=dim(2)
               if dx(i,j+1)~=0
       zjn(i,j)=boiln(i,j).*toyn(i,j)+bwatn(i,j).*twyn(i,j);
           a(counter,jjn-w(i,j+1))=zjn(i,j);
               end
           end
       zij(i,j)=-(boiln(i,j).*cop(i,j)+bwatn(i,j).*cwp(i,j)+...
           zin(i,j)+zip(i,j)+zjp(i,j)+zjn(i,j));
       
     if poi(i,j+1)<=0
        ip=i;
    else
        ip=i-1;
    end
    if poi(i+2,j+1)<=0
        in=i;
    else
        in=i+1;
    end
    if poi(i+1,j)<=0
        jp=j;
    else
        jp=j-1;
    end
    if poi(i+1,j+2)<=0
        jn=j;
    else
        jn=j+1;
    end
    qaqu;
           %%%%%%%%%%%%%%%%%%matrix a and b%%%%%%%%%%%%%%
dens=dens./144;
a(counter,position)=zij(i,j);              
b(counter,1)=-(boiln(i,j).*cop(i,j)+bwatn(i,j).*cwp(i,j)).*poil(i,j)-...
             (boiln(i,j).*qo(i,j)+bwatn(i,j).*qw(i,j)+qaqu(i,j))+bwatn(i,j).*(...
twxp(i,j).*(pcow(ip,j)-pcow(i,j))+...
twxn(i,j).*(pcow(in,j)-pcow(i,j))+...
twyp(i,j).*(pcow(i,jp)-pcow(i,j))+...
twyn(i,j).*(pcow(i,jn)-pcow(i,j)))+...
((boil(i,j).*toxp(i,j)).*dens(1)+(bwat(i,j).*twxp(i,j)).*dens(2)).*(z(ip,j)-z(i,j))+...
((boil(i,j).*toxn(i,j)).*dens(1)+(bwat(i,j).*twxn(i,j)).*dens(2)).*(z(in,j)-z(i,j))+...
((boil(i,j).*toyp(i,j)).*dens(1)+(bwat(i,j).*twyp(i,j)).*dens(2)).*(z(i,jp)-z(i,j))+...
((boil(i,j).*toyn(i,j)).*dens(1)+(bwat(i,j).*twyn(i,j)).*dens(2)).*(z(i,jn)-z(i,j));        
   end
end 
ppoilnup=inv(a)*b;
o=1;
for j=1:dim(2)
   for i=1:dim(1)
        if dx(i,j)==0
           o=o+1;
           continue
        end 
    position=i+(j-1)*dim(1)-o+1 ;
    poilnup(i,j)=ppoilnup(position,1);
    poilnup;
    err(i,j)=poiln(i,j)-poilnup(i,j);
    erro=err(i,j)+erro;
    error=erro./(dim(1).*dim(2));
   end
end
poiln=poilnup;
error;
                   end%while
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%saturation calculation%%%%%%%%%%%%%%%%%%
satwn=zeros(dim(1),dim(2));
saton=zeros(dim(1),dim(2));
for i=1:dim(1)
    for j=1:dim(2)
    if poi(i,j+1)<=0
        ip=i;
    else
        ip=i-1;
    end
    if poi(i+2,j+1)<=0
        in=i;
    else
        in=i+1;
    end
    if poi(i+1,j)<=0
        jp=j;
    else
        jp=j-1;
    end
    if poi(i+1,j+2)<=0
        jn=j;
    else
        jn=j+1;
    end
    if cww(i,j)==0
        satwn(i,j)=satw(i,j);
         if satwn(i,j)>1
           satw(i,j)=1;
       end
    else
        satwn(i,j)=satw(i,j)+(1./cww(i,j)).*((...
        twxp(i,j).*((poiln(ip,j)-poiln(i,j))-(pcow(ip,j)-pcow(i,j))-(dens(2).*(z(ip,j)-z(i,j))))+...
        twxn(i,j).*((poiln(in,j)-poiln(i,j))-(pcow(in,j)-pcow(i,j))-(dens(2).*(z(in,j)-z(i,j))))+...
        twyp(i,j).*((poiln(i,jp)-poiln(i,j))-(pcow(i,jp)-pcow(i,j))-(dens(2).*(z(i,jp)-z(i,j))))+...
        twyn(i,j).*((poiln(i,jn)-poiln(i,j))-(pcow(i,jn)-pcow(i,j))-(dens(2).*(z(i,jn)-z(i,j)))))-...
        cwp(i,j).*(poiln(i,j)-poil(i,j))+qw(i,j)+qaqu(i,j)./bwatn(i,j));
       if satwn(i,j)>=1
           satwn(i,j)=1;
       end
        saton(i,j)=1-satwn(i,j);
        if dx(i,j)==0
        satwn(i,j)=0;saton(i,j)=0;
        end
    end
    end
end
%%%%%%%%%%%%%pcow calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcown(i,j)=0;
for i=1:dim(1)
    for j=1:dim(2)
        if satwn(i,j)>pc(sizepc(1),1)
            pcown(i,j)==0;
        else
        pcown(i,j)=interp1(pc(:,1),pc(:,2),satwn(i,j),'linear'); 
        end
        pwatn(i,j)=poiln(i,j)-pcown(i,j);
        if dx(i,j)==0
           pwatn(i,j)=0;pcown(i,j)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%material balance check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imbo=0;
imbw=0;
eimbcw=0;
eimbco=0;
qostep=0;
qwstep=0;
von=0;
qaqustep=0;
for i=1:dim(1)
    for j=1:dim(2)
        if dx(i,j)==0
            imbo=0;
            imbw=0;
            eimbcw=0;
            eimbco=0;
            continue
        else
            if or(poi(i,j+1)==-1,poi(i+2,j+1)==-1)
               qaqu(i,j)=1.127.*1e-3.*(7000-pwatn(i,j)).*landaw(i,j).*((2.*kx(i,j).*dy(i,j).*dz(i,j))./(dx(i,j)));
            end
            if or(poi(i+1,j)==-1,poi(i+1,j+2)==-1)
               qaqu(i,j)=1.127.*1e-3.*(7000-pwatn(i,j)).*landaw(i,j).*((2.*ky(i,j).*dx(i,j).*dz(i,j))./(dy(i,j)))+qaqu(i,j);
            end
            if pwatn(i,j)>7000
               qaqu(i,j)=0;
            end
        qostep=qo(i,j)+qostep;
        qwstep=qw(i,j)+qwstep;
        qaqustep=qaqu(i,j)+qaqustep;
        von=(((dx(i,j).*dy(i,j).*dz(i,j))./5.615).*(por(i,j).*(sato(i,j)./boil(i,j))))+von;        
         if sato(i,j)==0
           continue
        end
        imbo=((((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*((por(i,j).*sato(i,j))./boil(i,j)))./...
                  ((((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*((porn(i,j).*saton(i,j))./boiln(i,j)))-qo(i,j)))+imbo;
        imbw=((((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*((por(i,j).*satw(i,j))./bwat(i,j)))./...
                  ((((dx(i,j).*dy(i,j).*dz(i,j))./(5.615.*time(1))).*((porn(i,j).*satwn(i,j))./bwatn(i,j)))-qw(i,j)-qaqu(i,j)./bwat(i,j)))+imbw;
        eimbco=(abs(imbo-1)).*100;
        eimbcw=(abs(imbw-1)).*100;
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%recovery calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qocum=qostep.*time(1)+qocum;
qwcum=qwstep.*time(1)+qwcum;
qaqucum=qaqustep.*time(1)+qaqucum;
rec=((voi-von)/voi)*100;
%%%%%%%%%%%%%%%%%%%%%output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qo
qw
qostep
qwstep
eimbco
eimbcw
pwatn
pcown
poiln
saton
satwn
qaqu
pwf
rec
qaqucum
qocum
qwcum
%%%%%%%%%%%%%%%%%%%%%%%%%update parameters%%%%%%%%%%%%%%%%%%%%%%%
satw=satwn;
sato=saton;
poil=poiln;
pwat=pwatn;
pcow=pcown;
    end%time

                          
                              
                  
                              
                              
                              
                              
