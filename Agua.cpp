#include <iostream>
#include "Agua.h"
#include <vector>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>

#define IX(i,j) ((i)+(N+2)*(j))
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

Agua::Agua(){

    N = 64;
    dt = 1.0f;
    diff = 0.0f;
    visc = 0.0f;
    force = 0.2f;
    source = 100.0f;
    allocate_data();
    dist_terreno = 1.5;
}

void Agua::desenharDensidade(){
    glBegin(GL_QUADS);
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        glColor3f(dens[IX(i+1,j+1)],dens[IX(i+1,j+1)],1.0/*dens[IX(i+1,j+1)]*/);
        //glNormal3f(0.0,-1.0,0.0);
        glVertex3f(i*float(10.0f/(N)) - 5.0f,dist_terreno,j*float(10.0f/(N)) - 5.0f);
        glVertex3f(i*float(10.0f/(N)) - 5.0f,dist_terreno,(j+1)*float(10.0f/(N)) - 5.0f);
        glVertex3f((i+1)*float(10.0f/(N)) - 5.0f,dist_terreno,(j+1)*float(10.0f/(N)) - 5.0f);
        glVertex3f((i+1)*float(10.0f/(N)) - 5.0f,dist_terreno,j*float(10.0f/(N)) - 5.0f);
    }
    /*glVertex3f(-5.0,3,-5.0);
    glVertex3f(5.0,3,-5.0);
    glVertex3f(5.0,3,5.0);
    glVertex3f(-5.0,3,5.0);*/
    glEnd();
}

void Agua::desenharDensidade2D(){
    float d00, d01, d10, d11;

    glBegin(GL_QUADS);
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
            d00 = dens[IX(i+1,j+1)];
            d01 = dens[IX(i+1,j+2)];
            d10 = dens[IX(i+2,j+1)];
            d11 = dens[IX(i+2,j+2)];

        glColor3f ( d00, d00, d00 );
        //glColor3f(dens[IX(i+1,j+1)],dens[IX(i+1,j+1)],dens[IX(i+1,j+1)]);
        glVertex2f(i*float(20.0f/(N)) - 10.0f,j*float(20.0f/(N)) - 10.0f);
         glColor3f ( d10, d10, d10 );
        glVertex2f(i*float(20.0f/(N)) - 10.0f,(j+1)*float(20.0f/(N)) - 10.0f);
        glColor3f ( d11, d11, d11 );
        glVertex2f((i+1)*float(20.0f/(N)) - 10.0f,(j+1)*float(20.0f/(N)) - 10.0f);
       glColor3f ( d01, d01, d01 );
        glVertex2f((i+1)*float(20.0f/(N)) - 10.0f,j*float(20.0f/(N)) - 10.0f);
    }
    /*glVertex3f(-5.0,3,-5.0);
    glVertex3f(5.0,3,-5.0);
    glVertex3f(5.0,3,5.0);
    glVertex3f(-5.0,3,5.0);*/
    glEnd();
}

void Agua::desenharVelocidade(){
    glBegin(GL_QUADS);
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        glColor3f(1.0f,1.0f,0.0f);
        glVertex3f(i*float(10.0f/(N)) - 5.0f,0.5,j*float(10.0f/(N)) - 5.0f);
        glVertex3f(i*float(10.0f/(N)) - 5.0f,0.5,(j+1)*float(10.0f/(N)) - 5.0f);
        glVertex3f((i+1)*float(10.0f/(N)) - 5.0f,0.5,(j+1)*float(10.0f/(N)) - 5.0f);
        glVertex3f((i+1)*float(10.0f/(N)) - 5.0f,0.5,j*float(10.0f/(N)) - 5.0f);
    }
    glEnd();
}

void Agua::add_source(int N, float * x, float * s, float dt){
   int i, size=(N+2)*(N+2);
   for ( i=0 ; i<size ; i++ )
   {
        x[i] += dt*s[i];
   }
}

void Agua::advect(int N, int b, float * d, float * d0, float * u, float * v, float dt){
   int i, j, i0, j0, i1, j1;
   float x, y, s0, t0, s1, t1, dt0;

   dt0 = dt*N;
   FOR_EACH_CELL
      x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
      if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
      if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
      s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
      d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
                s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
   END_FOR
   set_bnd ( N, b, d );
}

void Agua::advectDens(int N, int b, vector<Vetor> &d, vector<Vetor> &d0, float * u, float * v, float dt){
   int i, j, i0, j0, i1, j1;
   float x, y, s0, t0, s1, t1, dt0;

   dt0 = dt*N;
   FOR_EACH_CELL
      x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
      if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
      if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
      s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
      d[IX(i,j)].y = s0*(t0*d0[IX(i0,j0)].y+t1*d0[IX(i0,j1)].y)+
                s1*(t0*d0[IX(i1,j0)].y+t1*d0[IX(i1,j1)].y);
   END_FOR
   set_bnd_Dens ( N, b, d );
}

void Agua::project(int N, float * u, float * v, float * p, float * div){
   int i, j;

   FOR_EACH_CELL
      div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
      p[IX(i,j)] = 0;
   END_FOR
   set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

   lin_solve ( N, 0, p, div, 1, 4 );

   FOR_EACH_CELL
      u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
      v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
   END_FOR
   set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

void Agua::diffuse(int N, int b, float * x, float * x0, float diff, float dt){
   float a=dt*diff*N*N;
   lin_solve ( N, b, x, x0, a, 1+4*a );
}

void Agua::diffuseDens(int N, int b, vector<Vetor> &x, vector<Vetor> &x0, float diff, float dt){
   float a=dt*diff*N*N;
   lin_solve_Dens ( N, b, x, x0, a, 1+4*a );
}

void Agua::lin_solve_Dens(int N, int b, vector<Vetor> &x, vector<Vetor> &x0, float a, float c){
   int i, j, k;

   for ( k=0 ; k<20 ; k++ ) {
      FOR_EACH_CELL
         x[IX(i,j)].y = (x0[IX(i,j)].y + a*(x[IX(i-1,j)].y+x[IX(i+1,j)].y+x[IX(i,j-1)].y+x[IX(i,j+1)].y))/c;
      END_FOR
      set_bnd_Dens ( N, b, x );
   }
}

void Agua::lin_solve(int N, int b, float * x, float * x0, float a, float c){
   int i, j, k;

   for ( k=0 ; k<20 ; k++ ) {
      FOR_EACH_CELL
         x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
      END_FOR
      set_bnd ( N, b, x );
   }
}

void Agua::set_bnd( int N, int b, float * x){
   int i;

   for ( i=1 ; i<=N ; i++ ) {
      x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
      x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
      x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
      x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
   }
   x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
   x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
   x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
   x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
}

void Agua::set_bnd_Dens( int N, int b, vector<Vetor> &x){
   int i;

   for ( i=1 ; i<=N ; i++ ) {
      x[IX(0  ,i)].y = b==1 ? -x[IX(1,i)].y : x[IX(1,i)].y;
      x[IX(N+1,i)].y = b==1 ? -x[IX(N,i)].y : x[IX(N,i)].y;
      x[IX(i,0  )].y = b==2 ? -x[IX(i,1)].y : x[IX(i,1)].y;
      x[IX(i,N+1)].y = b==2 ? -x[IX(i,N)].y : x[IX(i,N)].y;
   }
   x[IX(0  ,0  )].y = 0.5f*(x[IX(1,0  )].y+x[IX(0  ,1)].y);
   x[IX(0  ,N+1)].y = 0.5f*(x[IX(1,N+1)].y+x[IX(0  ,N)].y);
   x[IX(N+1,0  )].y = 0.5f*(x[IX(N,0  )].y+x[IX(N+1,1)].y);
   x[IX(N+1,N+1)].y = 0.5f*(x[IX(N,N+1)].y+x[IX(N+1,N)].y);
}

void Agua::vel_step(int N, float * u, float * v, float * u0, float * v0, float visc, float dt){
   add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
   SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
   SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
   project ( N, u, v, u0, v0 );
   SWAP ( u0, u ); SWAP ( v0, v );
   advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
   project ( N, u, v, u0, v0 );
   SWAP ( u0, u ); SWAP ( v0, v );
}

void Agua::dens_step(int N, vector<Vetor> &x, vector<Vetor> &x0, float * u, float * v, float diff, float dt){
   //add_source ( N, x, x0, dt );
   //setDens(100,30,30);
   x.swap(x0); diffuseDens ( N, 0, x, x0, diff, dt );
   x.swap(x0); advectDens ( N, 0, x, x0, u, v, dt );
}

void Agua::roda_agua(Terreno *terreno){

   vel_step ( N, u, v, u_prev, v_prev, visc, dt );
   dens_step ( N, terreno->alturas, terreno->alturas0, u, v, diff, dt );

   if(isnan(dens[IX(30,30)]))
   {
      /*cout<<dens[IX(30,30)]<<endl;
      cout<<dens_prev[IX(30,30)]<<endl;
      cout<<u[IX(30,30)]<<endl;
      cout<<v[IX(30,30)]<<endl;
      cout<<u_prev[IX(30,30)]<<endl;
      cout<<v_prev[IX(30,30)]<<endl;*/

   }
}

void Agua::free_data ( void )
{
   if ( u ) free ( u );
   if ( v ) free ( v );
   if ( u_prev ) free ( u_prev );
   if ( v_prev ) free ( v_prev );
   if ( dens ) free ( dens );
   if ( dens_prev ) free ( dens_prev );
}

void Agua::clear_data ( void )
{
   int i, size=(N+2)*(N+2);

   for ( i=0 ; i<size ; i++ ) {
      u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
   }
}

int Agua::allocate_data ( void )
{
   int size = (N+2)*(N+2);

   u        = (float *) malloc ( size*sizeof(float) );
   v        = (float *) malloc ( size*sizeof(float) );
   u_prev      = (float *) malloc ( size*sizeof(float) );
   v_prev      = (float *) malloc ( size*sizeof(float) );
   dens     = (float *) malloc ( size*sizeof(float) );
   dens_prev   = (float *) malloc ( size*sizeof(float) );

   if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev ) {
      fprintf ( stderr, "cannot allocate data\n" );
      return ( 0 );
   }

   return ( 1 );
}

void Agua::setDens(float value, int i, int j){
    dens[IX(i,j)] = value;
}
void Agua::setU(float value, int i, int j){
    u_prev[IX(i,j)] = value;
}
void Agua::setV(float value, int i, int j){
    v_prev[IX(i,j)] = value;
}

Agua::~Agua(){
    clear_data();
    free_data();
}
