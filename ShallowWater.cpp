#include "ShallowWater.h"
#include <iostream>
#include <time.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#include <algorithm>

#define IX(i,j) ((i)+(N)*(j))

ShallowWater::ShallowWater(int N, Terreno *terreno, float dt)
{

    this->N = N;
    this->terreno = terreno;
    alturaAgua = vector<Vetor>(N*N);
    alturaAguaVelha = vector<Vetor>(N*N);
    flow = vector<Flow>(N*N);
    u = vector<float>(N*N);
    v = vector<float>(N*N);

    Ix = 10.0/N;
    Iy = 10.0/N;
    inicializaAlturas();
    inicializaVelocidadeAgua();
}

void ShallowWater::inicializaAlturas()
{
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
        {
            alturaAgua[IX(i,j)].x = float(i*(10.0f/float(N))) - 5.0f;
            alturaAgua[IX(i,j)].y = 0;//rand()%3;
            alturaAgua[IX(i,j)].z = float(j*(10.0f/float(N))) - 5.0f;
            alturaAguaVelha[IX(i,j)].x = alturaAgua[IX(i,j)].x;
            alturaAguaVelha[IX(i,j)].y = 0;
            alturaAguaVelha[IX(i,j)].z = alturaAgua[IX(i,j)].z;
            flow[IX(i,j)].L = 0;
            flow[IX(i,j)].R = 0;
            flow[IX(i,j)].T = 0;
            flow[IX(i,j)].B = 0;
            flow[IX(i,j)].K = 0;
        }
}

void ShallowWater::atualizaAlturaVelha(){
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        alturaAguaVelha[IX(i,j)].y = alturaAgua[IX(i,j)].y;
        /*flow[IX(i,j)].L = 1;
        flow[IX(i,j)].R = 1;
        flow[IX(i,j)].T = 1;
        flow[IX(i,j)].B = 1;
        flow[IX(i,j)].K = 1;*/
    }
}

void ShallowWater::atualizaAlturas(int i, int j){
    alturaAgua[IX(i,j)].y = alturaAgua[IX(i,j)].y + terreno->alturas[IX(i,j)].y;
}

int ShallowWater::inicializaVelocidadeAgua(){
   //int size = (N)*(N);

   //u        = (float *) malloc ( size*sizeof(float) );
   //v        = (float *) malloc ( size*sizeof(float) );

   for(int i=0; i<N; i++)
   for(int j=0; j<N; j++)
   {
       u[IX(i,j)] = 0;
       v[IX(i,j)] = 0;
   }


   return 1;
}


void ShallowWater::free_data ( void )
{
//   if ( u ) free ( u );
 //  if ( v ) free ( v );
}


void ShallowWater::desenhar()
{

    glBegin(GL_TRIANGLES);

    float alturaDesenhada;
    float dim = 10.0f/float(N);

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        if(alturaAgua[IX(i,j)].y>0)
        {
                alturaDesenhada = alturaAgua[IX(i,j)].y + terreno->alturas[IX(i,j)].y;
                glColor3f(alturaAgua[IX(i,j)].y,alturaAgua[IX(i,j)].y,0.5);
                glNormal3f(0,1,0);
                glVertex3f(alturaAgua[IX(i,j)].x-dim/2,alturaDesenhada,alturaAgua[IX(i,j)].z+dim/2);
                glVertex3f(alturaAgua[IX(i,j)].x-dim/2,alturaDesenhada,alturaAgua[IX(i,j)].z-dim/2);
                glVertex3f(alturaAgua[IX(i,j)].x+dim/2,alturaDesenhada,alturaAgua[IX(i,j)].z+dim/2);

                glColor3f(alturaAgua[IX(i,j)].y,alturaAgua[IX(i,j)].y,0.5);
                glNormal3f(0,1,0);
                glVertex3f(alturaAgua[IX(i,j)].x+dim/2,alturaDesenhada,alturaAgua[IX(i,j)].z+dim/2);
                glVertex3f(alturaAgua[IX(i,j)].x-dim/2,alturaDesenhada,alturaAgua[IX(i,j)].z-dim/2);
                glVertex3f(alturaAgua[IX(i,j)].x+dim/2,alturaDesenhada,alturaAgua[IX(i,j)].z-dim/2);
        }
    }

    glEnd();
}

void ShallowWater::waterIncrement(int i,int j,float r)
{
    alturaAgua[IX(i,j)].y += r*dt*kr;
}

void ShallowWater::rainDrop(float r)
{
    //srand ( time(NULL) );
    int i = rand()%N;
    int j = rand()%N;
    waterIncrement(i,j,r);
    //atualizaAlturas(i,j);
}

void ShallowWater::fixedFont(int i, int j, float r)
{
    waterIncrement(i,j,r);
    //atualizaAlturas(i,j);
}


void ShallowWater::waterEvaporation(){

    float ke = 0.15;

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        alturaAgua[IX(i,j)].y = max(0.0f, alturaAgua[IX(i,j)].y*(1-ke*dt));

        /*if(alturaAgua[IX(i,j)].y>0.05f){
            alturaAgua[IX(i,j)].y = 0.0f;
        }*/
    }
}

void ShallowWater::evaporateAllWater(){
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        alturaAgua[IX(i,j)].y = 0.0f;
        alturaAguaVelha[IX(i,j)].y = 0.0f;
        /*if(alturaAgua[IX(i,j)].y>0.05f){
            alturaAgua[IX(i,j)].y = 0.0f;
        }*/
    }
}

void ShallowWater::inflowFlux(int i, int j)
{

}

void ShallowWater::outflowFlux(int i, int j)
{
    float diferenca;
    float rho=1.0f;
    float l=1.0f;
    float A =2.0f;
    float g=9.81f;
    float fluxFactor = rho*dt*A*g/l;
    float alturaAtual = alturaAgua[IX(i,j)].y + terreno->alturas[IX(i,j)].y;
    if(i>=1)
    {
        diferenca = alturaAtual - (alturaAgua[IX(i-1,j)].y + terreno->alturas[IX(i-1,j)].y);
        flow[IX(i,j)].L = max(0.0f, flow[IX(i,j)].L + diferenca*fluxFactor);
    }
    else
    flow[IX(i,j)].L = 0.0f;


    if(i<=N-2)
    {
        diferenca = alturaAtual - (alturaAgua[IX(i+1,j)].y + terreno->alturas[IX(i+1,j)].y);
        flow[IX(i,j)].R = max(0.0f, flow[IX(i,j)].R + diferenca*fluxFactor);
    }
    else
    flow[IX(i,j)].R = 0.0f;


    if(j>=1)
    {
        diferenca = alturaAtual - (alturaAgua[IX(i,j-1)].y + terreno->alturas[IX(i,j-1)].y);
        flow[IX(i,j)].B = max(0.0f, flow[IX(i,j)].B + diferenca*fluxFactor);
    }
    else
    flow[IX(i,j)].B = 0.0f;


    if(j<=N-2)
    {
        diferenca = alturaAtual - (alturaAgua[IX(i,j+1)].y + terreno->alturas[IX(i,j+1)].y);
        flow[IX(i,j)].T = max(0.0f, flow[IX(i,j)].T + diferenca*fluxFactor);
    }
    else
    flow[IX(i,j)].T = 0.0f;

    float somaFluxos = flow[IX(i,j)].L + flow[IX(i,j)].R + flow[IX(i,j)].T + flow[IX(i,j)].B;

    if(somaFluxos<=0)
        somaFluxos = 1;

    flow[IX(i,j)].K = min(1.0f, (alturaAgua[IX(i,j)].y*Ix*Iy)/(dt*somaFluxos));

    flow[IX(i,j)].L *= flow[IX(i,j)].K;
    flow[IX(i,j)].R *= flow[IX(i,j)].K;
    flow[IX(i,j)].B *= flow[IX(i,j)].K;
    flow[IX(i,j)].T *= flow[IX(i,j)].K;
}

void ShallowWater::newVolume(int i,int j){

    float inflow = (flow[IX(i+1,j)].L+flow[IX(i-1,j)].R+flow[IX(i,j+1)].B+flow[IX(i,j-1)].T);

    float outflow = (flow[IX(i,j)].L+flow[IX(i,j)].R+flow[IX(i,j)].B+flow[IX(i,j)].T);

    float dV = dt*(inflow - outflow);

    float alturaVelha = alturaAgua[IX(i,j)].y;
    alturaAgua[IX(i,j)].y += dV/(Ix*Iy);
    float alturaNova = alturaAgua[IX(i,j)].y;
    float d = 0.5f*(alturaNova + alturaVelha);

    if(d==0)d=1;

    u[IX(i,j)] = 0.5*(flow[IX(i-1,j)].R - flow[IX(i,j)].L + flow[IX(i,j)].R - flow[IX(i+1,j)].L)/(d*Iy);
    v[IX(i,j)] = 0.5*(flow[IX(i,j-1)].T - flow[IX(i,j)].B + flow[IX(i,j)].T - flow[IX(i,j+1)].B)/(d*Ix);
}

void ShallowWater::trataBordas(){

    for(int i=0; i<N; i++)
    {
        alturaAgua[IX(i,0)].y = alturaAgua[IX(i,1)].y;
        alturaAgua[IX(i,N-1)].y = alturaAgua[IX(i,N-2)].y;
        alturaAgua[IX(0,i)].y = alturaAgua[IX(1,i)].y;
        alturaAgua[IX(N-1,i)].y = alturaAgua[IX(N-2,i)].y;
    }


    for(int i=0; i<alturaAgua.size(); i++)
    {
        if(alturaAgua[i].y < 0)
           alturaAgua[i].y = 0;

        if(alturaAgua[i].y > 20)
           alturaAgua[i].y = 20;
    }

    alturaAgua[IX(0,0)].y = 0;
}

void ShallowWater::shallowWaterStep()
{

    /*for(int i=0; i<N-1; i++)
    for(int j=0; j<N-1; j++)
    {
        outflowFlux(i,j);
    }
    for(int i=0; i<N-1; i++)
    for(int j=0; j<N-1; j++)
    {
        newVolume(i,j);
    }/**/

    for(int i=0; i<N-1; i++)
    for(int j=0; j<N-1; j++)
    {
        outflowFlux(i,j);
        newVolume(i,j);

    }/**/
    trataBordas();

    atualizaAlturaVelha();
}

ShallowWater::~ShallowWater()
{
    free_data();
}
