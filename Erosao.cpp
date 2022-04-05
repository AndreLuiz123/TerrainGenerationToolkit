#include "Erosao.h"
#include <iostream>
#include "Agua.h"
#include <vector>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#define IX(i,j) ((i)+(N+2)*(j))
#define IXT(i,j) ((i)+(N)*(j))
#define PI 3.14159265

Erosao::Erosao(int N, Terreno *terreno, ShallowWater *agua, float kc, float kd, float ks){

    this->N=N;

    this->terreno = terreno;
    this->agua = agua;
    inicializaSedimentoSuspenso();

    this->dt = agua->dt;

    this->kc = kc;
    this->kd = kd;
    this->ks = ks;
    this->kt = 0.3;
    this->ka = 0.5;
}

void Erosao::inicializaSedimentoSuspenso(){
     //int size = (N)*(N);
     //sedimento_suspenso = (float *) malloc ( size*sizeof(float) );
     //sedimento_suspenso0 = (float *) malloc ( size*sizeof(float) );

    /* for(int i=0; i<N; i++)
     for(int j=0; j<N; j++)
     {
         sedimento_suspenso0[IXT(i,j)] = 0;
         sedimento_suspenso[IXT(i,j)] = 0;
     }*/
     sedimento_suspenso0 = vector<float>(N*N,0);
     sedimento_suspenso = vector<float>(N*N,0);

    /* for(int i=N/2-N/8; i<N/2+N/8; i++)
     for(int j=N/2-N/8; j<N/2+N/8; j++)
     {
         sedimento_suspenso0[IXT(i,j)] = 50;
         sedimento_suspenso[IXT(i,j)] = 50;
     }*/
}

normalize(Vetor *vec){

    float tamAtual = sqrt(vec->x*vec->x + vec->y*vec->y + vec->z*vec->z);

    vec->x = vec->x/tamAtual;
    vec->y = vec->y/tamAtual;
    vec->z = vec->z/tamAtual;
}

float dot(Vetor vec1, Vetor vec2){

    float prodEsc = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    float prodModulos = sqrt(vec1.x*vec1.x + vec1.y*vec1.y + vec1.z*vec1.z)*sqrt(vec2.x*vec2.x + vec2.y*vec2.y + vec2.z*vec2.z);

    return prodEsc/prodModulos;
}

void Erosao::roda_erosao(bool evapo){

    //cout<<"Deposita Dissolve: "<<endl<<sedimento_suspenso[IXT(N/2,N/2)]<<" "<<agua->alturaAgua[IXT(N/2,N/2)].y<<endl;
    //cout<<"Transporta Sedimento: "<<endl<<sedimento_suspenso0[IXT(N/2+N/8+1,N/2)]<<" "<<agua->u[IXT(N/2,N/2)]<<endl;
    agua->shallowWaterStep();
    deposita_dissolve();
    sedimento_suspenso.swap(sedimento_suspenso0);
    transportaSedimento();
    guardaGradeSedimentoSuspenso();
    sedimento_suspenso.swap(sedimento_suspenso0);
    thermalErosion();
    if(evapo) agua->waterEvaporation();
}

void Erosao::roda_erosao_hidraulica(bool evapo){
    agua->shallowWaterStep();
    deposita_dissolve();
    sedimento_suspenso.swap(sedimento_suspenso0);
    transportaSedimento();
    guardaGradeSedimentoSuspenso();
    sedimento_suspenso.swap(sedimento_suspenso0);
    //guardaGradeSedimentoSuspenso();
    if(evapo)
    agua->waterEvaporation();
}

void Erosao::aeolianErosion(){
    vento->roda_agua(terreno);
}

void Erosao::roda_erosao_termica(){
    thermalErosion();
}

void Erosao::desenha_sedimento(){

    float dim = 10.0f/float(N);
    float alturaDesenhada=0.01f;
    glBegin(GL_QUADS);
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        if(sedimento_suspenso[IXT(i,j)]>1000)
        {
            //cout<<"IX("<<i<<","<<j<<"): "<<sedimento_suspenso[IX(i,j)]<<endl;
            break;
        }
        if(sedimento_suspenso[IXT(i,j)]>0)
        {
            alturaDesenhada = terreno->alturas[IXT(i,j)].y + 0.5;
            glColor3f(0.5,sedimento_suspenso[IXT(i,j)],sedimento_suspenso[IXT(i,j)]);
            glNormal3f(0,1,0);
            glVertex3f(float(i*(10.0f/float(N))) - 5.0f-dim/2,alturaDesenhada,float(j*(10.0f/float(N))) - 5.0f+dim/2);
            glVertex3f(float(i*(10.0f/float(N))) - 5.0f-dim/2,alturaDesenhada,float(j*(10.0f/float(N))) - 5.0f-dim/2);
            glVertex3f(float(i*(10.0f/float(N))) - 5.0f+dim/2,alturaDesenhada,float(j*(10.0f/float(N))) - 5.0f-dim/2);
            glVertex3f(float(i*(10.0f/float(N))) - 5.0f+dim/2,alturaDesenhada,float(j*(10.0f/float(N))) - 5.0f+dim/2);
        }

    }
    glEnd();
}

void Erosao::deposita_dissolve(){

    //cout<<"Deposita dissolve"<<sedimento_suspenso[IXT(N/2,N/2)]<<" "<<agua->alturaAgua[IXT(N/2,N/2)].y<<endl;
    float uV;
    float vV;

    for(int i=1; i<N; i++)
    for(int j=1; j<N; j++)
    {
       uV = agua->u[IXT(i,j)];
       vV = agua->v[IXT(i,j)];

       float b[4];
        b[0] = terreno->alturas[IXT(i,j)].y - terreno->alturas[IXT(i+1,j)].y;
        b[1] = terreno->alturas[IXT(i,j)].y - terreno->alturas[IXT(i-1,j)].y;
        b[2] = terreno->alturas[IXT(i,j)].y - terreno->alturas[IXT(i,j+1)].y;
        b[3] = terreno->alturas[IXT(i,j)].y - terreno->alturas[IXT(i,j-1)].y;

        float localTiltY = max(max(b[0],b[1]),max(b[2],b[3]));
        float localTiltX = (10.0/float(N));
        float sinAlpha = localTiltY/(sqrt(localTiltY*localTiltY + localTiltX*localTiltX));/**/

        /*Vetor normal;
        //if(i-1>0 || i+1<N-1)
        normal.x = terreno->alturas[IXT(i+1,j)].y - terreno->alturas[IXT(i-1,j)].y;

        //float xY = terreno->alturas[IXT(i+1,j)].y - terreno->alturas[IXT(i-1,j)].y;
        //float zY = terreno->alturas[IXT(i,j+1)].y - terreno->alturas[IXT(i,j-1)].y;

        normal.y = 1;
        //if(j-1>0 || j+1<N-1)
        normal.z = terreno->alturas[IXT(i,j+1)].y - terreno->alturas[IXT(i,j-1)].y;
        normalize(&normal);
        Vetor up;
        up.x=0;
        up.y=1;
        up.z=0;
        float cosa = dot(normal,up);/**/
        ///float sinAlpha = sqrt(1 - cosa*cosa);//sin(acos(cosa));/**/
        sinAlpha = std::max(sinAlpha,0.1f);
        //sinAlpha = std::max(sinAlpha,0.1f);


        float C= kc*sinAlpha*sqrtf((uV*uV)+(vV*vV));

        if(agua->alturaAgua[IXT(i,j)].y>0)
        {
            if(C<=sedimento_suspenso[IXT(i,j)]){
                depositaSedimento(i,j,sedimento_suspenso[IXT(i,j)] - C);
            }else {
                dissolveSedimento(i,j,C - sedimento_suspenso[IXT(i,j)]);
            }
        }else{
            depositaSedimento(i,j,sedimento_suspenso[IXT(i,j)]);
        }
        }
}

void Erosao::dissolveSedimento(int i,int j, float ds){

    //cout<<ds<<endl;
    if(terreno->alturas[IXT(i,j)].y>0)
    {
        sedimento_suspenso[IXT(i,j)] += ks*ds*0.001;
        terreno->alturas[IXT(i,j)].y -= ks*ds*0.001;
    }

    if(terreno->alturas[IXT(i,j)].y<0)
        terreno->alturas[IXT(i,j)].y = 0;

    //if(ds<2.5)
        //cout<<"ks: "<<ks<<" ds: "<<ds<<" ds*ks"<<ds*ks<<" sedimento suspenso: "<<sedimento_suspenso[IXT(i,j)]<<endl;
}

void Erosao::depositaSedimento(int i,int j, float ds){

    if(sedimento_suspenso[IXT(i,j)]>0)
    {
        sedimento_suspenso[IXT(i,j)] -= kd*ds*0.001;
        terreno->alturas[IXT(i,j)].y += kd*ds*0.001;
        //cout<<"Deposita ds:"<<ds<<" kd:"<<kd<<" kd*ds:"<<ds*kd<<endl;
    }


    if(sedimento_suspenso[IXT(i,j)]<0)
        sedimento_suspenso[IXT(i,j)] = 0;
}


void Erosao::transportaSedimento(){

   int i, j, i0, j0, i1, j1;
   float x, y, s0, t0, s1, t1;

    float dt0 = dt*N;

    for(int i=1; i<N-1; i++)
    for(int j=1; j<N-1; j++)
    {


        //x = float(i) - dt*agua->u[IXT(i,j)];
        //y = float(j) - dt*agua->v[IXT(i,j)];

        x = agua->alturaAgua[IXT(i,j)].x - dt0*agua->u[IXT(i,j)];
        y = agua->alturaAgua[IXT(i,j)].z - dt0*agua->v[IXT(i,j)];

      if (x<-4.5f) x=-4.5f; if (x>4.5f) x=4.5f;
      if (y<-4.5f) y=-4.5f; if (y>4.5f) y=4.5f;

      float n = N;
      float gs = 10.0f/n;

      float x2 = x + 5;
      float y2 = y + 5;

      i0 = int(x2/gs);
      i1 = i0+1;

      j0 = int(y2/gs);
      j1 = j0+1;

      s1 = (x2-i0*gs)/(gs); s0 = 1-s1; t1 = (y2-j0*gs)/(gs); t0 = 1-t1;

      sedimento_suspenso[IXT(i,j)] = (s0*(t0*sedimento_suspenso0[IXT(i0,j0)]+t1*sedimento_suspenso0[IXT(i0,j1)])+
                s1*(t0*sedimento_suspenso0[IXT(i1,j0)]+t1*sedimento_suspenso0[IXT(i1,j1)]));
    //else
    //sedimento_suspenso[IXT(i,j)] = 0;
        //if(t1>0)
            //cout<<"s0:"<<s0<<" s1:"<<s1<<" t0:"<<t0<<" t1:"<<t1<<" sed susp 0:"<<sedimento_suspenso0[IXT(i1,j0)]<<endl;

    }

    for ( i=0 ; i<N; i++ ) {
      sedimento_suspenso[IXT(0  ,i)] = sedimento_suspenso[IXT(1,i)];
      sedimento_suspenso[IXT(N-1,i)] = sedimento_suspenso[IXT(N-2,i)];
      sedimento_suspenso[IXT(i,0  )] = sedimento_suspenso[IXT(i,1)];
      sedimento_suspenso[IXT(i,N-1)] = sedimento_suspenso[IXT(i,N-2)];
   }

   sedimento_suspenso[IXT(0  ,0  )] = 0.5f*(sedimento_suspenso[IXT(1,0  )]+sedimento_suspenso[IXT(0  ,1)]);
   sedimento_suspenso[IXT(0  ,N-1)] = 0.5f*(sedimento_suspenso[IXT(1,N-1)]+sedimento_suspenso[IXT(0  ,N-2)]);
   sedimento_suspenso[IXT(N-1,0  )] = 0.5f*(sedimento_suspenso[IXT(N-2,0  )]+sedimento_suspenso[IXT(N-1,1)]);
   sedimento_suspenso[IXT(N-1,N-1)] = 0.5f*(sedimento_suspenso[IXT(N-2,N-1)]+sedimento_suspenso[IXT(N-1,N-2)]);
}

void Erosao::thermalErosion(){

    terreno->alturas0 = terreno->alturas;
    float d = (10.0/float(N));
    float a = d*d;

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        float b[8];
        float newdS[8];


        b[0] = i<N-1?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i+1,j)].y:-1.0f;
        b[1] = i<N-1 && j<N-1?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i+1,j+1)].y:-1.0f;
        b[2] = j<N-1?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i,j+1)].y:-1.0f;
        b[3] = i>0 && j<N-1?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i-1,j+1)].y:-1.0f;
        b[4] = i>0?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i-1,j)].y:-1.0f;
        b[5] = i>0 && j>0?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i-1,j-1)].y:-1.0f;
        b[6] = j>0?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i,j-1)].y:-1.0f;
        b[7] = i<N-1 && j>0?terreno->alturas0[IXT(i,j)].y - terreno->alturas0[IXT(i+1,j-1)].y:-1.0f;

        float H = max(max(max(max(b[0],b[1]),b[2]),b[3]),max(max(max(b[4],b[5]),b[6]),b[7]));

        float sum;
        for(int k=0; k<8; k++)
            b[k]>0?sum+=b[k]:sum+=0;

        float dS = kt*dt*H*0.5f*a;

        for(int k=0; k<8; k++)
        {
            float talus = b[k]/d;

            if(b[k]>0 && talus>1)//(b[k]/d))
            {
                //cout<<"erodiu ("<<i<<","<<j<<")"<<endl;
                newdS[k] = dS*b[k]/sum;
                terreno->alturas[IXT(i,j)].y -=newdS[k];
            }else
                newdS[k] = 0;
        }

        //terreno->alturas[IXT(i,j)].y -= dS;
        terreno->alturas[IXT(i+1,j)].y+=newdS[0];
        terreno->alturas[IXT(i+1,j+1)].y+=newdS[1];
        terreno->alturas[IXT(i,j+1)].y+=newdS[2];
        terreno->alturas[IXT(i-1,j+1)].y+=newdS[3];
        terreno->alturas[IXT(i-1,j)].y+=newdS[4];
        terreno->alturas[IXT(i-1,j-1)].y+=newdS[5];
        terreno->alturas[IXT(i,j-1)].y+=newdS[6];
        terreno->alturas[IXT(i+1,j-1)].y+=newdS[7];
    }
}

void Erosao::transportaSedimento2(){

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        float sedimento = 0.0f;
        int den = 0;
        //sedimento_suspenso[IXT(i,j)]=0;
        if(i>0)
        {
            //sedimento+=sedimento_suspenso0[IXT(i-1,j)];
            if(sedimento_suspenso0[IXT(i,j)] < sedimento_suspenso0[IXT(i-1,j)] && agua->alturaAgua[IXT(i-1,j)].y>0)
            {
                sedimento_suspenso[IXT(i,j)]+=sedimento_suspenso0[IXT(i-1,j)]/4;
                sedimento_suspenso[IXT(i-1,j)]-=sedimento_suspenso0[IXT(i-1,j)]/4;
            }

            den++;
        }
        if(i+1<N)
        {
            //sedimento+=sedimento_suspenso0[IXT(i-1,j)];
            if(sedimento_suspenso0[IXT(i,j)] < sedimento_suspenso0[IXT(i+1,j)]&& agua->alturaAgua[IXT(i+1,j)].y>0)
            {
                sedimento_suspenso[IXT(i,j)]+=sedimento_suspenso0[IXT(i+1,j)]/4;
                sedimento_suspenso[IXT(i+1,j)]-=sedimento_suspenso0[IXT(i+1,j)]/4;
            }
            //sedimento+=sedimento_suspenso0[IXT(i+1,j)];
            den++;
        }
        if(j>0)
        {
            //sedimento+=sedimento_suspenso0[IXT(i-1,j)];
            if(sedimento_suspenso0[IXT(i,j)] < sedimento_suspenso0[IXT(i,j-1)]&& agua->alturaAgua[IXT(i,j-1)].y>0)
            {
                sedimento_suspenso[IXT(i,j)]+=sedimento_suspenso0[IXT(i,j-1)]/4;
                sedimento_suspenso[IXT(i,j-1)]-=sedimento_suspenso0[IXT(i,j-1)]/4;
            }
            //sedimento+=sedimento_suspenso0[IXT(i,j-1)];
            den++;
        }
        if(j+1<N)
        {
            //sedimento+=sedimento_suspenso0[IXT(i-1,j)];
            if(sedimento_suspenso0[IXT(i,j)] < sedimento_suspenso0[IXT(i,j+1)]&& agua->alturaAgua[IXT(i,j+1)].y>0)
            {
                sedimento_suspenso[IXT(i,j)]+=sedimento_suspenso0[IXT(i,j+1)]/4;
                sedimento_suspenso[IXT(i,j+1)]-=sedimento_suspenso0[IXT(i,j+1)]/4;
            }
            //sedimento+=sedimento_suspenso0[IXT(i,j+1)];
            den++;
        }
    }
}

void Erosao::guardaGradeSedimentoSuspenso(){
    //for(int i=0; i<sedimento_suspenso0.size(); i++)
        sedimento_suspenso0 = sedimento_suspenso;

}

Erosao::~Erosao(){
    //terreno->~terreno();
    //terreno = NULL;
    //delete terreno;
    //agua->~Agua();
    //agua = NULL;
    //delete agua;
    //if ( sedimento_suspenso ) free ( sedimento_suspenso );
    //if ( sedimento_suspenso0 ) free ( sedimento_suspenso0 );
}
