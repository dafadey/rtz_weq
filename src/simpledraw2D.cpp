#ifdef _WIN32
#  define WINDOWS_LEAN_AND_MEAN
#  include <windows.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <pthread.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
//#include <rendercheck_gl.h>

//using namespace std;
static GLuint* texture;
static int*    texUPD;
static int*    texSZ;
static float*    texMAG;
static float*    texPOS;
static float** texturesDataArrs;
static int* close_win;
static pthread_mutex_t				close_win_mutex;

static int XSize = 512;
static int YSize = 512;
static int XnewSize = XSize;
static int YnewSize = YSize;
static int Nx;
static int Ny;
static int count;
static int side_max;
static int xmouse, ymouse;
static float pdx_drag, pdy_drag;
static int inited=0;

static float red(float val, float max)
{
	if(val>0)	return sqrt(val/max);
	else return (val/max)*(val/max)*0.35;
}
static float green(float val, float max)
{
	if(val<0) return (val/max)*(val/max)*0.75+sqrt(-val/max)*0.25;
	else	return (val/max)*(val/max);
}
static float blue(float val, float max)
{
	if(val<0)	return sqrt(-val/max);
	else return 0.0;
}
static void mouseXY(int button, int state, int x, int y)
{
	xmouse=x;
	ymouse=y;
	pdx_drag=0;
	pdy_drag=0;
	int wheel=0;
//	if((button==GLUT_LEFT_BUTTON)&(state==GLUT_DOWN))
//		printf("GL: x=%f, y=%f\n",float(x)/float(XnewSize),float(y)/float(YnewSize));
	if ((button==3)&(state==GLUT_DOWN)) wheel=+1;
	if ((button==4)&(state==GLUT_DOWN)) wheel=-1;
//	if((button==GLUT_RIGHT_BUTTON)&(state==GLUT_DOWN))	printf("GL: XSize=%d\n",XSize);
	
	if(wheel!=0)
	{
		float xa=float(x)/float(XnewSize);
		float ya=float(y)/float(YnewSize);
		int ii=floor(xa*float(side_max));
		int jj=floor(ya*float(side_max));
		int number=ii+jj*side_max;
		float pmag=texMAG[number];
		texMAG[number]+=float(wheel)/13.0;
		float lxa=xa*float(side_max)-float(ii);
		float lya=ya*float(side_max)-float(jj);
		texPOS[number*2]+=lxa*(pmag-texMAG[number]);
		texPOS[number*2+1]+=lya*(pmag-texMAG[number]);
	}
}


static void mouse_drag(int x, int y)
{

	float xa=float(xmouse)/float(XnewSize);
	float ya=float(ymouse)/float(YnewSize);
	int ii=floor(xa*float(side_max));
	int jj=floor(ya*float(side_max));
	int number=ii+jj*side_max;

	int dx=x-xmouse;
	int dy=y-ymouse;
//	printf("GL: drag x=%d y=%d\n",dx,dy);
	float dxa=float(dx)/float(XnewSize)*float(side_max)*texMAG[number];
	float dya=float(dy)/float(YnewSize)*float(side_max)*texMAG[number];
	texPOS[number*2]-=dxa-pdx_drag;
	texPOS[number*2+1]-=dya-pdy_drag;	
	pdx_drag=dxa;
	pdy_drag=dya;
}

static void win_sz(int x, int y)
{
		glViewport(0,0,x,y);
		XnewSize=x;
		YnewSize=y;
//		printf("GL (Resize): XSize=%d\n",XSize);
}

static void 
display(void)
{
	glutSwapBuffers();
  glMatrixMode (GL_MODELVIEW);
  glDisable(GL_DEPTH_TEST);
  glClearColor(.3, .3, .3, 0);
  glClear(GL_COLOR_BUFFER_BIT);
	int ii,jj;
	int x=0;
	int y=0;
	int dx=floor(float(XSize)/float(side_max));
	int dy=floor(float(YSize)/float(side_max));

	for(jj=0;jj<side_max;jj++)
	{
		x=0;
		for(ii=0;ii<side_max;ii++)
		{
			int number=ii+jj*side_max;
			if(number<count)
			{
				float mag=texMAG[number];
				float xtex=texPOS[number*2];
				float ytex=texPOS[number*2+1];
				glBindTexture(GL_TEXTURE_2D, texture[number]);
				glBegin(GL_QUADS);
				glTexCoord2f(xtex, ytex);
				glVertex2f(x+1, y+1);
				glTexCoord2f(xtex+1.0*mag, ytex);
				glVertex2f(x+dx-1, y+1);
				glTexCoord2f(xtex+1.0*mag, ytex+1.0*mag);
				glVertex2f(x+dx-1, y+dy-1);
				glTexCoord2f(xtex, ytex+1.0*mag);
				glVertex2f(x+1, y+dy-1);
				glEnd();
			}
			x+=dx;
		}
		y+=dy;
	}
}

static void idle()
{
 	  int i;
    for(i=0;i<count;i++)
    {
        if(texUPD[i]==0)
        {
		    glBindTexture(GL_TEXTURE_2D, texture[i]);
#ifdef GL_VERSION_1_1
		    glTexImage2D(GL_TEXTURE_2D, 0, 4,texSZ[i*2],texSZ[i*2+1], 0, GL_RGBA, GL_FLOAT, texturesDataArrs[i]);
#else
		    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB,texSZ[i*2],texSZ[i*2+1], 0, GL_RGBA, GL_FLOAT, texturesDataArrs[i]);
#endif
		    texUPD[i]=1;
		}
        
    }
	glutPostRedisplay();
	if(close_win[0]==1)
	{
		delete[] texture;
		delete[] texUPD;
		delete[] texSZ;
		delete[] texMAG;
		delete[] texPOS;
		for(i=0;i<count;i++) delete[] texturesDataArrs[i];
		delete[] texturesDataArrs;
		delete[] close_win;
		exit(0);
	}
//	pthread_mutex_unlock(&close_win_mutex);
}

static pthread_t					GLloop_th_id;
void* GLloop(void*)
{

//	printf("GLloop: activated\n");
	int i,j,k;
	texturesDataArrs=new float*[count];
	texUPD=new int[count];
	texSZ=new int[count*2];
	texMAG=new float[count];
	texPOS=new float[count*2];
	for(k=0;k<count;k++)
	{
		texSZ[k*2]=Nx;
		texSZ[k*2+1]=Ny;
		texMAG[k]=1.0;
		texPOS[k*2]=0.0;
		texPOS[k*2+1]=0.0;
	}
	for(i=0;i<count;i++) texUPD[i]=0;
	for(i=0;i<count;i++) texturesDataArrs[i]=new float[texSZ[i*2]*texSZ[i*2+1]*4];
	texture=new GLuint[count];

	for(k=0;k<count;k++)
	{
		for(j=0;j<texSZ[k*2+1];j++)
		{
		  for(i=0;i<texSZ[k*2];i++)
		  {
			 	float x=float(i-Nx/2)/float(Nx);
				float y=float(j-Ny/2)/float(Ny);
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4]=exp(-(x*x+y*y)*13.0*(1+k));
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4+1]=0.0;
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4+2]=0.0;
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4+3]=0.0;
		  }
		}
		//texturesDataArrs[k][((64+25)*128+64+33)*4+1]=1.0;
		//texturesDataArrs[k][((64-33)*128+64+47)*4+2]=1.0;
	}
	char* argv;//='0';
	int argc=1;
	glutInit(&argc, &argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(XSize, YSize);
//	int main_w=
	glutCreateWindow("FADEY");
  
	glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho (0, XSize, YSize, 0, -1, 1);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);
	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glutReshapeFunc(win_sz);
//	glutSpecialFunc(processSpecialKeys); 

	glutMouseFunc(mouseXY);
	glutMotionFunc(mouse_drag);
//	glewInit();

	int ii, jj;
	side_max=-floor(-sqrt(float(count)));
//  glActiveTexture(GL_TEXTURE0);
  glGenTextures(count, texture);
	for(jj=0;jj<side_max;jj++)
	{
		for(ii=0;ii<side_max;ii++)
		{
			int number=ii+jj*side_max;
			if(number<count)
			{
				//printf("making texture number %d\n",number);
				glBindTexture(GL_TEXTURE_2D, texture[number]);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
#ifdef GL_VERSION_1_1
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, texSZ[number*2], texSZ[number*2+1], 0, GL_RGBA, GL_FLOAT, texturesDataArrs[number]);
#else
				glTexImage2D(GL_TEXTURE_2D, 0, 4, texSZ[number*2], texSZ[number*2+1], 0, GL_RGBA, GL_FLOAT, texturesDataArrs[number]);
#endif
				texUPD[number]=1;				
			}
		}
	}
	//printf("side_max=%d\n",side_max);
	pthread_mutex_unlock(&close_win_mutex);
	glutMainLoop();
//	printf("fadey_draw: bye\n");
	return 0;
}

void fadey_init(int Nx_, int Ny_, int count_)
{
	if(inited==0)
	{
		inited=1;
		Nx=Nx_;
		Ny=Ny_;
		count=count_;
		pthread_create(&GLloop_th_id, NULL, &GLloop, NULL);
		close_win=new int[1];
		close_win[0]=0;
		pthread_mutex_init(&close_win_mutex,NULL);
		pthread_mutex_lock(&close_win_mutex);
		pthread_mutex_lock(&close_win_mutex);
	//	printf("init: bye\n");
	}
}


void fadey_draw(float* DataArr, int Nx_, int Ny_, int count_)
{
	if(inited==1)
	{  
		int i,j;
		int k=count_;
		if ((Nx_!=texSZ[k*2])||(Ny_!=texSZ[k*2+1]))
		{
			delete[] texturesDataArrs[k];
			texSZ[k*2]=Nx_;
			texSZ[k*2+1]=Ny_;
			texturesDataArrs[k]=new float[texSZ[k*2]*texSZ[k*2+1]*4];
//			texturesDataArrs[k][(j*texSZ[k*2]+i)*4]=DataArr[j*texSZ[k*2]+i];
			texturesDataArrs[k]=DataArr;
	//		printf("texture %d new size is (%d x %d)\n",k,texSZ[k*2],texSZ[k*2+1]);
		}
		float max=0.0;
		float min=0.0;
		for(j=0;j<texSZ[k*2+1];j++)
		{
			for(i=0;i<texSZ[k*2];i++)
			{
				if(DataArr[j*texSZ[k*2]+i]>max) max=DataArr[j*texSZ[k*2]+i];
				if(DataArr[j*texSZ[k*2]+i]<max) min=DataArr[j*texSZ[k*2]+i];
			}
		}
		if(-min>max) max=-min;
		for(j=0;j<texSZ[k*2+1];j++)
		{
			for(i=0;i<texSZ[k*2];i++)
			{
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4]=red(DataArr[j*texSZ[k*2]+i],max);
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4+1]=green(DataArr[j*texSZ[k*2]+i],max);
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4+2]=blue(DataArr[j*texSZ[k*2]+i],max);
				texturesDataArrs[k][(j*texSZ[k*2]+i)*4+3]=0.0;
			}
		}
		texUPD[k]=0;
  }  
}

void fadey_close()
{
	if(inited==1)
	{
		close_win[0]=1;
		pthread_mutex_lock(&close_win_mutex);
	}
}

