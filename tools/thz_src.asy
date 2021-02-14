settings.outformat="pdf";
import graph;
import libimaging_debug;
defaultpen(linewidth(0.1));

struct anotations	{
	int filenum;
	string anotation;
	void operator init(int i,string s){
		this.filenum=i;
		this.anotation=s;
	}
}

//file fin=input("../result.dat");
file fin=input(flname);
string a[]=fin;

bool work=true;
picture p;
int cl;
while(cl<a.length)
{
	//read to new array
	int N=(int)a[cl];
	cl+=1;
	string samp=a[cl];
	cl+=1;
	string sphase=a[cl];
	cl+=1;
	string scurv=a[cl];
	cl+=1;
	real[] s=new real[N]; 
	for(int i=0;i<N;++i)
	{
		s[i]=(real)a[cl];
		cl+=1;
	}
	//
	guide g;
	real amp=(real)samp;
	real phase=(real)sphase;
	real curv=(real)scurv;
	//find max
	real max=0.0;
	int imax;
	for(int i=0;i<N;++i)
	{
		if(max<s[i])
		{
			max=s[i];
			imax=i;
		}
	}
	//find width
	real wp=max/2;
	int iwbeg=0;
	int iwend=0;
	for(int i=0;i<N-1;++i)
	{
		if((s[i]<wp)&(s[i+1]>wp)) iwbeg=i; 
		if((s[i]>wp)&(s[i+1]<wp)) iwend=i; 
	}
	dot(p,(imax,s[imax]),red);
	dot(p,(iwbeg,s[iwbeg]),blue);
	dot(p,(iwend,s[iwend]),blue);
	for(int i=0;i<N;++i)	g=g--(i,s[i]);
	draw(p,g,phase/1.57*red);
	
	anotations[] toplot={anotations(imax,"Maximum"),anotations(iwbeg,"left"),anotations(iwend,"right")};
	//draw field on axis for given extra points
	//test_w2amp=0.1_w2pahse=0.0_w2curv=1.0
	for(int i=0;i<toplot.length;++i)
	{
		string filename="../test_w2amp="+samp+"_w2pahse="+sphase+"_w2curv="+scurv+"/pics/res/field"+string(toplot[i].filenum)+".dat";
		write(filename+" "+toplot[i].anotation);
		file fsigin=input(filename);
		real sig[]=fsigin;

		filename="../test_w2amp="+samp+"_w2pahse="+sphase+"_w2curv="+scurv+"/pics/res/plasma"+string(toplot[i].filenum)+".dat";
		write(filename+" "+toplot[i].anotation);
		file fthzin=input(filename);
		real THZ[]=fthzin;


		//read
		int Nsig=1024;
		int Nspec=floor(Nsig/2);
		real[] sigw=new real[Nsig];
		real[] sig2w=new real[Nsig];
		real[] spec=new real[Nspec];
		//(s>0.000001)?pow(s*s+0.0000000000000001,-0.3125)*exp(-1.0/(fabs(s)+0.000000000000000000001)):0.0;
		real[] n=new real[Nsig];
		real[] thzs=new real[Nsig];
		real[] THZs=new real[Nsig]; // from file
		{
			real np=0.0;
			for(int j=0;j<Nsig;++j)
			{
				real s=sig[j];
				sigw[j]=s;
				sig2w[j]=s;
				np+=(abs(s)>0.000001)?(s*s+0.0000000000000001)^(-0.3125)*exp(-1.0/(abs(s)+0.000000000000000000001)):0.0;
				n[j]=np;
				thzs[j]=np*s;
				THZs[j]=THZ[j]; //reading source from file
			}
		}
		//filter
		int ffreq=70;
		int lowpassfreq=10;
		qqFFT_freal(Nsig,sigw);
		qqFFT_freal(Nsig,sig2w);
		qqFFT_freal(Nsig,thzs);
		for(int j=0;j<Nspec;++j)
		{
			spec[j]=sqrt(sigw[j*2]*sigw[j*2]+sigw[j*2+1]*sigw[j*2+1]);
			sigw[j*2]=(j<ffreq)?sigw[j*2]:0.0;
			sigw[j*2+1]=(j<ffreq)?sigw[j*2+1]:0.0;
			sig2w[j*2]=(j>ffreq)?sig2w[j*2]:0.0;
			sig2w[j*2+1]=(j>ffreq)?sig2w[j*2+1]:0.0;
			thzs[j*2]=exp(-j^2/lowpassfreq^2)*thzs[j*2];
			thzs[j*2+1]=exp(-j^2/lowpassfreq^2)*thzs[j*2+1];
		}
		qqFFT_freal_1(Nsig,sigw);
		qqFFT_freal_1(Nsig,sig2w);
		qqFFT_freal_1(Nsig,thzs);
		//normalize n
		{
			real max=0.0;
			for(int j=0;j<Nsig;++j) max=(max>sigw[j])?max:sigw[j];
			for(int j=0;j<Nsig;++j) n[j]=n[j]/n[Nsig-1]*max;
			real thzsmax=0.0;
			for(int j=0;j<Nsig;++j) thzsmax=(thzsmax>thzs[j])?thzsmax:thzs[j];
			for(int j=0;j<Nsig;++j) thzs[j]=thzs[j]/thzsmax*max;
			real THZsmax=0.0;
			for(int j=0;j<Nsig;++j) THZsmax=(THZsmax>THZs[j])?THZsmax:THZs[j];
			for(int j=0;j<Nsig;++j) THZs[j]=THZs[j]/THZsmax*max;
		}
		//draw
		guide gsigw,gsig2w,gspec,gn,gthzs,gTHZs;
		picture subp,subps;
		for(int j=0;j<Nsig;++j)
		{
			gsigw=gsigw--(j,sigw[j]);
			gsig2w=gsig2w--(j,sig2w[j]);
			gn=gn--(j,n[j]);
			gthzs=gthzs--(j,thzs[j]);
			gTHZs=gTHZs--(j,THZs[j]);
		}
		for(int j=0;j<Nspec;++j) gspec=gspec--(j,spec[j]);
		draw(subp,gsig2w,blue);
		draw(subp,gsigw,red);
		draw(subp,gn);
		draw(subp,gthzs,0.43*green);
		draw(subp,gTHZs,0.43*green+dashed);
		draw(subp,(0,0)--(Nsig-1,0));
		draw(subps,gspec);
		dot(subps,(ffreq,0.0));
		dot(subps,(lowpassfreq,0.0));
		real X=6cm;
		real Y=1.5cm;
		xaxis(subp,"$\tau$",BottomTop,LeftTicks("$\tiny{%.4g}$"),true);
		yaxis(subp,"$E_{opt}$",LeftRight,RightTicks("$\tiny{%.4g}$"),true);
		size(subp,X,Y,point(subp,SW),point(subp,NE));
		label(subp,toplot[i].anotation,point(subp,NE),SW);
		add(subp.fit(),(i*(X+2cm),(7cm+Y+2cm)));

		xaxis(subps,"$\omega$",BottomTop,LeftTicks("$\tiny{%.4g}$"),true);
		yaxis(subps,"$E_{opt}$",LeftRight,RightTicks("$\tiny{%.4g}$"),true);
		size(subps,X,Y,point(subps,SW),point(subps,NE));
		add(subps.fit(),(i*(X+2cm),(7cm+2*Y+4cm)));
	}
}
xaxis(p,"$z$\,(arb. units)",BottomTop,LeftTicks("$\tiny{%.4g}$"),true);
yaxis(p,"$S_{THz}$",LeftRight,RightTicks("$\tiny{%.4g}$"),true);
size(p,21cm,7cm,point(p,SW),point(p,NE));
add(p.fit(),(0,0));
