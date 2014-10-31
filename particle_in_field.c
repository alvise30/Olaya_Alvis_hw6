#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>

#define CARGA_E 1.602176E-19 // Carga del electron
#define MASA_E 1.672622E-27  // Masa del electron
#define VEL_C 299792458.0    // Velocidad de la luz
#define PI 3.141593
#define B_0 3.0E-5          
#define RT 6378100.0
#define TMIN 0.0
#define TMAX 100.0

double bx_prima(double posx, double posy, double posz, double velx, double vely, double velz, double gamma);
double by_prima(double posx, double posy, double posz, double velx, double vely, double velz, double gamma);
double bz_prima(double posx, double posy, double posz, double velx, double vely, double velz, double gamma);


int main(int argc, char **argv){
  double h;
  double x_old, y_old, z_old;
  double xprime_old, yprime_old, zprime_old;
  double x_new, y_new, z_new;
  double xprime_new, yprime_new, zprime_new;
  double x1, y1,z1, t1;
  double x2, y2,z2, t2;
  double x3, y3,z3, t3;
  double xprime1, yprime1,zprime1;
  double xprime2, yprime2,zprime2;
  double xprime3, yprime3,zprime3;
  double kx1, ky1, kz1;
  double kx2, ky2, kz2;
  double kx3, ky3, kz3;
  double kx4, ky4, kz4;
  double kx5, ky5, kz5;
  double E0 ,E0j,alpha,v0, gamma,v1;
  int i,n_puntos;

//constantes
E0 = atof(argv[1]);
alpha = atof(argv[2]);
E0j = E0*CARGA_E*1000000;
gamma =1+E0j/(MASA_E*VEL_C*VEL_C);
v0=VEL_C*sqrt(1-(1/pow(gamma,2)));
h=(2*PI*gamma*MASA_E/(CARGA_E*B_0))/50;

n_puntos=(TMAX-TMIN)/h;

//Crear archivo
FILE *fileout;
char filename[100000];
 sprintf(filename, "trayectoria_%d_%d.dat", (int)E0, (int)alpha);
fileout = fopen(filename,"w");

//condiciones iniciales
x_old=2*RT;
y_old=0.0;
z_old=0.0;
xprime_old=0;
yprime_old=v0*sin(alpha* PI / 180.0);
zprime_old=v0*cos(alpha * PI / 180.0);
fprintf(fileout, "%f\t%f\t%f\t%f\t%f\t%f\t%f \n",x_old,y_old,z_old,xprime_old,yprime_old,zprime_old,v0);

//rungekutteira
for(i=1;i<n_puntos;i++){
  kx1=bx_prima(x_old,y_old,z_old,xprime_old,yprime_old,zprime_old,gamma);
  ky1=by_prima(x_old,y_old,z_old,xprime_old,yprime_old,zprime_old,gamma);
  kz1=bz_prima(x_old,y_old,z_old,xprime_old,yprime_old,zprime_old,gamma);
//primer paso
  x1=x_old+h*xprime_old/2.0;
  y1=y_old+h*yprime_old/2.0;
  z1=z_old+h*zprime_old/2.0;
  xprime1=xprime_old+h*kx1/2.0;
  yprime1=yprime_old+h*ky1/2.0;
  zprime1=zprime_old+h*kz1/2.0;
  kx2=bx_prima(x1,y1,z1,xprime1,yprime1,zprime1,gamma);
  ky2=by_prima(x1,y1,z1,xprime1,yprime1,zprime1,gamma);
  kz2=bz_prima(x1,y1,z1,xprime1,yprime1,zprime1,gamma);
//segundo paso
  x2=x_old+h*xprime1/2.0;
  y2=y_old+h*yprime1/2.0;
  z2=z_old+h*zprime1/2.0;
  xprime2=xprime_old+h*kx2/2.0;
  yprime2=yprime_old+h*ky2/2.0;
  zprime2=zprime_old+h*kz2/2.0;
  kx3=bx_prima(x2,y2,z2,xprime2,yprime2,zprime2,gamma);
  ky3=by_prima(x2,y2,z2,xprime2,yprime2,zprime2,gamma);
  kz3=bz_prima(x2,y2,z2,xprime2,yprime2,zprime2,gamma);
//tercer paso
  x3=x_old+h*xprime2;
  y3=y_old+h*yprime2;
  z3=z_old+h*zprime2;
  xprime3=xprime_old+h*kx3;
  yprime3=yprime_old+h*ky3;
  zprime3=zprime_old+h*kz3;
  kx4=bx_prima(x3,y3,z3,xprime3,yprime3,zprime3,gamma);
  ky4=by_prima(x3,y3,z3,xprime3,yprime3,zprime3,gamma);
  kz4=bz_prima(x3,y3,z3,xprime3,yprime3,zprime3,gamma);
//cuarto paso
  kx5=(1.0/6.0)*(kx1 + 2.0*kx2 + 2.0*kx3 + kx4);
  ky5=(1.0/6.0)*(ky1 + 2.0*ky2 + 2.0*ky3 + ky4);
  kz5=(1.0/6.0)*(kz1 + 2.0*kz2 + 2.0*kz3 + kz4);
  xprime_new=xprime_old+h*kx5;
  yprime_new=yprime_old+h*ky5;
  zprime_new=zprime_old+h*kz5;
  x_new=x_old+h*(xprime_new+xprime_old)/2;
  y_new=y_old+h*(yprime_new+yprime_old)/2;
  z_new=z_old+h*(zprime_new+zprime_old)/2;
  v1=sqrt(pow(xprime_new,2)+pow(yprime_new,2)+pow(zprime_new,2));
  xprime_old=xprime_new;
  yprime_old=yprime_new;
  zprime_old=zprime_new;
  x_old=x_new;
  y_old=y_new;
  z_old=z_new;
fprintf(fileout, "%f\t%f\t%f\t%f\t%f\t%f\t%f \n",x_new,y_new,z_new,xprime_new,yprime_new,zprime_new,v1);
}
return 0;
}

double bx_prima(double posx, double posy, double posz, double velx, double vely, double velz, double gamma){
  double r;
  r=pow(posx*posx+posy*posy+posz*posz,0.5);
  return CARGA_E*(-B_0*pow(RT,3.0))*(vely*(2.0*posz*posz-posx*posx-posy*posy)-velz*(3.0*posy*posz))/(MASA_E*gamma*pow(r,5.0));
  }

double by_prima(double posx, double posy, double posz, double velx, double vely, double velz, double gamma){
  double r;
  r=pow(posx*posx+posy*posy+posz*posz,0.5);
  return CARGA_E*(-B_0*pow(RT,3.0))*(velz*(3.0*posx*posz)-velx*(2.0*posz*posz-posx*posx-posy*posy))/(MASA_E*gamma*pow(r,5.0));
 }

double bz_prima(double posx, double posy, double posz, double velx, double vely, double velz, double gamma){
  double r;
  r=pow(posx*posx+posy*posy+posz*posz,0.5);
  return CARGA_E*(-B_0*pow(RT,3.0))*(velx*(3.0*posy*posz)-vely*(3.0*posx*posz))/(MASA_E*gamma*pow(r,5.0));  
  }
