settings.outformat="pdf";
import graph;
defaultpen(linewidth(0.3));

file fin=input("rp_sharp30w.dat");
real a[]=fin;
int Nw=30;
int Nth=floor(a.length/Nw);
write("number of theta points is ",Nth);
picture p;
for(int i=1;i<Nw;i+=1)
{
	real MAX=0.0;
	for(int j=0;j<Nth;++j) MAX=(MAX>a[j*Nw+i])?MAX:a[j*Nw+i];
	guide rp;
	for(int j=0;j<Nth;++j) rp=rp--(j/Nth*0.1,a[j*Nw+i]/MAX);
	draw(p,rp,(i+1)/Nw*red);
	write("freq="+string(i*1.17)+" THz");
	i=i+floor(i/2);
}
xaxis(p,"$\theta$ [rad.]",BottomTop,LeftTicks("$\tiny{%.4g}$"),true);
yaxis(p,"$W_{THz}$",LeftRight,RightTicks("$\tiny{%.4g}$"),true);
size(p,7cm,5cm,point(p,SW),point(p,NE));
add(p.fit());


