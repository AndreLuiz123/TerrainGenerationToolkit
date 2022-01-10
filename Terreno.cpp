#include "Terreno.h"

#include <iostream>
#include <time.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#include <fstream>
#include "Erosao.h"
#define IX(i,j) ((i)+(N)*(j))
#define PI 3.14159265

Terreno::Terreno(){

    this->N = 64;
    alturas = vector<Vetor>(N*N);
    inicializaAlturas();
    hillAlgorithm(5);

}

Terreno::Terreno(int N, int numberHills){
    srand ( time(NULL) );
    this->N = N;
    alturas = vector<Vetor>(N*N);
    alturasY = vector<float>(N*N);
    alturas0 = vector<Vetor>(N*N);
    alturasY0 = vector<float>(N*N);
    inicializaAlturas();
    if(numberHills>0)
    //faulting(numberHills);
    //hillAlgorithm(numberHills);
    diamond(numberHills);
}

Vetor calculaNormal(Vetor p1, Vetor p2, Vetor p3){
    Vetor u,v;
    Vetor normal;

    u.x = p2.x - p1.x;
    u.y = p2.y - p1.y;
    u.z = p2.z - p1.z;

    v.x = p3.x - p1.x;
    v.y = p3.y - p1.y;
    v.z = p3.z - p1.z;

    normal.x = u.y*v.z - v.y*u.z;
    normal.y = u.z*v.x - v.z*u.x;
    normal.z = u.x*v.y - v.x*u.y;

    return normal;
}

normalize2(Vetor *vec){

    float tamAtual = sqrt(vec->x*vec->x + vec->y*vec->y + vec->z*vec->z);

    vec->x = vec->x/tamAtual;
    vec->y = vec->y/tamAtual;
    vec->z = vec->z/tamAtual;
}

void Terreno::desenhar(){

Vetor normal;
glBegin(GL_TRIANGLES);
    for(int i=0; i<N-1; i++)
    for(int j=0; j<N-1; j++)
    {
        float r=0.0f, g=0.5f, b=0.0f;
        if(alturas[IX(i,j)].y>0.5)
        {
            r = 0.5f;
            g = 0.3f;
            b = 0.1f;
        }
        if(alturas[IX(i,j)].y>3.5)
        {
            r = 0.8f;
            g = 0.8f;
            b = 0.8f;
        }
        /*
        r = alturas[IX(i,j)].y/5;
        g = 0.8f;
        b = alturas[IX(i,j)].y/10*/;

        normal = calculaNormal(alturas[IX(i,j)],alturas[IX(i,j+1)],alturas[IX(i+1,j)]);
        normalize2(&normal);
        glNormal3f(normal.x,normal.y,normal.z);
        glColor3f(r,g,b);
        glVertex3f(alturas[IX(i,j)].x,alturas[IX(i,j)].y,alturas[IX(i,j)].z);
        glVertex3f(alturas[IX(i,j+1)].x,alturas[IX(i,j+1)].y,alturas[IX(i,j+1)].z);
        glVertex3f(alturas[IX(i+1,j)].x,alturas[IX(i+1,j)].y,alturas[IX(i+1,j)].z);

        normal = calculaNormal(alturas[IX(i+1,j)],alturas[IX(i,j+1)],alturas[IX(i+1,j+1)]);
        normalize2(&normal);
        glNormal3f(normal.x,normal.y,normal.z);
        glColor3f(r,g,b);
        glVertex3f(alturas[IX(i+1,j)].x,alturas[IX(i+1,j)].y,alturas[IX(i+1,j)].z);
        glVertex3f(alturas[IX(i,j+1)].x,alturas[IX(i,j+1)].y,alturas[IX(i,j+1)].z);
        glVertex3f(alturas[IX(i+1,j+1)].x,alturas[IX(i+1,j+1)].y,alturas[IX(i+1,j+1)].z);

    }
glEnd();
glBegin(GL_POINTS);
    glColor3f(0.0,1.0,0.0);
glEnd();
}

void Terreno::inicializaAlturas(){
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        alturas[IX(i,j)].x = float(i*(10.0f/float(N))) - 5.0f;
        alturas[IX(i,j)].y = 0;//rand()%3;
        alturas[IX(i,j)].z = float(j*(10.0f/float(N))) - 5.0f;
        alturas0[IX(i,j)].x = float(i*(10.0f/float(N))) - 5.0f;
        alturas0[IX(i,j)].y = 0;//rand()%3;
        alturas0[IX(i,j)].z = float(j*(10.0f/float(N))) - 5.0f;
        alturasY[IX(i,j)] = 0;
        alturasY0[IX(i,j)] = 0;//rand()%3;
    }
}

float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

void Terreno::hillAlgorithm(int numberHills){

    int hillRadiusMin = N/8;
    int hillRadiusMax = N/5.33;


    float hillMinHeight = 0.0;
    float hillMaxHeight = 59.0;

    float maxValueGenerated=-1;
    float minValueGenerated=-1;

    //srand ( time(NULL) );
    for(int h=0; h<numberHills; h++)
    {
        int centerCollum = rand()%N;
        //cout<<centerCollum<<endl;
        int centerRow = rand()%N;
        //cout<<centerRow<<endl;
        float hillHeight = RandomFloat(hillMinHeight,hillMaxHeight);
        int hillRadius = rand()%(hillRadiusMax - hillRadiusMin) + hillRadiusMin;

        for(int c = centerCollum - hillRadius; c<centerCollum + hillRadius; c++)
        for(int l = centerRow- hillRadius; l<centerRow + hillRadius; l++)
        {

            if(c<0 || c>=N || l<0 || l>=N)
                continue;

                float r2 = hillRadius*hillRadius;
                float x2x1 = centerCollum - c;
                float y2y1 = centerRow - l;
                float altura = r2 - x2x1*x2x1 - y2y1*y2y1;

                if(altura<0)
                    continue;

                float fator = altura/r2;

                alturas[IX(c,l)].y += hillHeight*fator;
                //if(alturas[IX(c,l)].y>20.0)
                    //alturas[IX(c,l)].y=20.0;

                if(maxValueGenerated==-1 || maxValueGenerated<alturas[IX(c,l)].y)
                    maxValueGenerated = alturas[IX(c,l)].y;

                if(minValueGenerated==-1 || minValueGenerated>alturas[IX(c,l)].y)
                    minValueGenerated = alturas[IX(c,l)].y;
        }
    }

    normalizeHeight(maxValueGenerated, minValueGenerated);
}

void Terreno::faultingDefaultFunction(int i, int j,float a, float b, float c, float displacement){

    if(a*alturas[IX(i,j)].x + b*alturas[IX(i,j)].z - c>0)
        alturas[IX(i,j)].y+=displacement;
    else
        alturas[IX(i,j)].y-=displacement;
}
void Terreno::faultingSinFunction(int i, int j,float a, float b, float c, float displacement){

    float dist = a*alturas[IX(i,j)].x + b*alturas[IX(i,j)].z - c;

        if(dist>0){
            if(dist<PI/2)
            alturas[IX(i,j)].y+=sin(dist)*displacement;
            else
            alturas[IX(i,j)].y+=displacement;
        }else{
            if(dist>-PI/2)
            alturas[IX(i,j)].y+=sin(dist)*displacement;
            else
            alturas[IX(i,j)].y-=displacement;
        }
}

void Terreno::faultingCosFunction(int i, int j,float a, float b, float c, float displacement){

    float dist = a*alturas[IX(i,j)].x + b*alturas[IX(i,j)].z - c;

    if(dist<PI/2 && dist>-PI/2)
        alturas[IX(i,j)].y+=cos(dist)*displacement;

}

void Terreno::faulting(int iterations){

    float w = alturas[IX(N-1,0)].x - alturas[IX(0,0)].x;
    float l = alturas[IX(0,N-1)].z - alturas[IX(0,0)].z;

    float maxValueGenerated=0;
    float minValueGenerated=0;
    for(int k=0; k<iterations; k++)
    {
        int v = rand()%360;
        float a = sin(v);
        float b = cos(v);
        float d = sqrt(w*w + l*l);
        float c = rand()%(int)d - d/2;
        float displacement = 0.3f;

        /*cout<<"c:"<<c<<endl;
        cout<<"v:"<<v<<endl;
        cout<<"a:"<<a<<endl;
        cout<<"b:"<<b<<endl;*/

        for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
        {
            faultingDefaultFunction(i,j,a,b,c,displacement);
            ///faultingSinFunction(i,j,a,b,c,displacement);
            ///faultingCosFunction(i,j,a,b,c,displacement);

            if(maxValueGenerated<alturas[IX(i,j)].y)
                maxValueGenerated = alturas[IX(i,j)].y;

            if(minValueGenerated>alturas[IX(i,j)].y)
                minValueGenerated = alturas[IX(i,j)].y;
        }
    }
    normalizeHeight(maxValueGenerated, minValueGenerated);
}

float Terreno::findMaxValue(bool maxOrMin){

    float maxValue=0;
    float minValue=alturas[0].y;

    for(int i=0; i<alturas.size(); i++)
    {
        if(alturas[i].y>maxValue)
            maxValue = alturas[i].y;
        if(alturas[i].y<minValue)
            minValue = alturas[i].y;
    }

    return maxOrMin?maxValue:minValue;
}

void Terreno::diamond_squareStep(int chunk, float dHeight, float r){
    cout<<"square"<<endl;
    for(int i=0; i<N-1; i+=chunk)
    for(int j=0; j<N-1; j+=chunk)
    {
        float randHeight =  RandomFloat(-dHeight/2, dHeight/2);
        //cout<<i<<" "<<j<<endl;
        alturas[IX(i+chunk/2,j+chunk/2)].y =
        (alturas[IX(i,j)].y + alturas[IX(i+chunk,j)].y +
        alturas[IX(i,j+chunk)].y + alturas[IX(i+chunk,j+chunk)].y)/4 + randHeight;
    }
}

void Terreno::diamond_diamondStep(int chunk, float dHeight, float r){
    cout<<"diamond"<<endl;
    int half = chunk/2;
    for(int i=0; i<N; i+=half)
    for(int j=(i+half)%chunk; j<N; j+=chunk)
    {
        float randHeight = RandomFloat(-dHeight/2, dHeight/2);
        float newHeight = 0;
        int cont=0;

        newHeight+=i>0?alturas[IX(i-half,j)].y:0;
        cont+=i>0?1:0;
        newHeight+=i<N-1?alturas[IX(i+half,j)].y:0;
        cont+=i<N-1?1:0;
        newHeight+=j>0?alturas[IX(i,j-half)].y:0;
        cont+=j>0?1:0;
        newHeight+=j<N-1?alturas[IX(i,j+half)].y:0;
        cont+=j<N-1?1:0;

        alturas[IX(i,j)].y = (newHeight/cont) + randHeight;
        //cout<<alturas[IX(i,j)]<<endl;
    }
}

void Terreno::diamond(int iterations){

    int chunk = (N/2);
    float r = chunk*0.01;
    float dHeight = 20;
    /*alturas[IX(0,0)].y = 40;//RandomFloat(0,40);///50+RandomFloat(-r,r);
    alturas[IX(N-1,0)].y = 40;//RandomFloat(0,40);///50+RandomFloat(-r,r);
    alturas[IX(0,N-1)].y = 40;//RandomFloat(0,40);///80+RandomFloat(-r,r);
    alturas[IX(N-1,N-1)].y = 40;//RandomFloat(0,40);///70+RandomFloat(r,r);*/
    alturas[IX((N-1)/2,(N-1)/2)].y = 20;
    //recursiveDiamond(iterations, 0, 0, N-1, 50);

    cout<<"aoba"<<endl;
    for(int i=0; i<iterations; i++)
    while(chunk>1){
        //\float half = chunck/2;
        diamond_squareStep(chunk,dHeight,r);
        diamond_diamondStep(chunk,dHeight,r);
        chunk/=2;
        dHeight/=2;
    }

    float maxVal = findMaxValue(true);
    float minVal = findMaxValue(false);

    normalizeHeight(maxVal,minVal);
}

void Terreno::recursiveDiamond(int iterations, int x0,int y0, int h, float dHeight){

    //cout<<iterations<<endl;
    if(iterations>0)
    {

        float randHeight =  RandomFloat(-dHeight/2, dHeight/2);
        float newHeight =
        (alturas[IX(x0,y0)].y + alturas[IX(x0,y0+h)].y + alturas[IX(x0+h,y0)].y + alturas[IX(x0+h,y0+h)].y)*0.25 + randHeight*0.5;

        randHeight =  RandomFloat(0.0f, dHeight);
        float newHeight1 = (alturas[IX(x0,y0)].y + alturas[IX(x0+h,y0)].y)*0.5 + dHeight*0.5;
        randHeight =  RandomFloat(0.0f, dHeight);
        float newHeight2 = (alturas[IX(x0,y0+h)].y + alturas[IX(x0+h,y0+h)].y)*0.5 + dHeight*0.5;
        randHeight =  RandomFloat(0.0f, dHeight);
        float newHeight3 = (alturas[IX(x0,y0)].y + alturas[IX(x0,y0+h)].y)*0.5 + dHeight*0.5;
        randHeight =  RandomFloat(0.0f, dHeight);
        float newHeight4 = (alturas[IX(x0+h,y0)].y + alturas[IX(x0+h,y0+h)].y)*0.5 + dHeight*0.5;


        alturas[IX(x0+h/2, y0+h/2)].y = newHeight;
        alturas[IX(x0+h/2, y0)].y = newHeight1;
        alturas[IX(x0+h/2, y0+h)].y = newHeight2;
        alturas[IX(x0, y0+h/2)].y = newHeight3;
        alturas[IX(x0+h, y0+h/2)].y = newHeight4;

        dHeight*=0.5;
        h/=2;
        recursiveDiamond(iterations-1,x0,y0,h,dHeight);
        recursiveDiamond(iterations-1,x0+h,y0,h,dHeight);
        recursiveDiamond(iterations-1,x0+h,y0+h,h,dHeight);
        recursiveDiamond(iterations-1,x0,y0+h,h,dHeight);
    }
}


void Terreno::perlinNoise(){

}

void Terreno::salvarTerreno(){
}

void Terreno::normalizeHeight(float maxH, float minH){

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        ///Normalize
        alturas[IX(i,j)].y = (alturas[IX(i,j)].y - minH)/(maxH - minH);
        ///Flatering
        alturas[IX(i,j)].y *= alturas[IX(i,j)].y;
        ///Put in the range 0 to 5
        alturas[IX(i,j)].y *= 5;

        //cout<<alturas[IX(i,j)].y<<endl;
    }
}

Terreno::~Terreno(){

}
