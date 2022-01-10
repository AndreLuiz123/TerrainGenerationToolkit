#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "Terreno.h"
#include "Agua.h"
#include "Erosao.h"
#include "ShallowWater.h"
#define IX(i,j) ((i)+(8)*(j))
int WIDTH=500, HEIGHT=500;

bool visual3D = true;
bool desenhaAgua3d = true;
bool rain = false;
bool panorama = true;

void display();
void reshape(int, int);
void timer(int);
void inicializaTerreno();
void keyboardFunction(unsigned char, int, int);
void pre_display_2D();
void pre_display_3D();


float x=-10.0,z =-10.0;
float velX=0.05, velZ=0.5;

float dt=0.008f;
int N=64;
Agua agua;
Terreno terreno(N, 122);
ShallowWater shallowWater(N, &terreno,dt);
Erosao erosao(N, &terreno, &shallowWater,1.0f, 0.4f,  0.7f);
string path = "C:\\Users\\Andre\\Documents\\TCC\\PIC C++\\ErosaoHidraulica3\\Arquivosalvo\\";

fstream file;

void saveData(Terreno *t){
    file.open(path+"terrenoSalvo.txt",ios::out);
    for(int i=0; i<t->alturas.size(); i++)
        file<<t->alturas[i].y<<endl;
}

void saveDataAndClose(Terreno *t){


    /*terrenoSalvo<<endl<<endl<<"Terreno ao final"<<endl<<endl;

    for(int i=0; i<t->alturas.size(); i++)
    {
        terrenoSalvo<<t->alturas[i].y<<endl;
    }
    //terrenoSalvo.close();*/
}
void init(){
    glClearColor(0.0,0.0,0.0,1.0);
    pre_display_3D();
    saveData(&terreno);
    //erosao.sedimento_suspenso[IX(N/2,N/2)] = 5;
    //erosao.sedimento_suspenso[IX(N/2,10)] = 0;
}


int main(int argc, char **argv){


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutInitWindowPosition(100,100);
    glutInitWindowSize(500,500);


    glutCreateWindow("Erosao Hidraulica");
    glutKeyboardFunc(keyboardFunction);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(1,timer,0);
    init();

    glutMainLoop();
}

void codigoErosaoHidraulica1(){

    if(visual3D)
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        glTranslatef(0,0,-16.0);
        glRotated(60.0,1.0,0.0,0.0);
        glRotated(x,0.0,1.0,0.0);

        x++;

        float matSpec[] ={1.0f,1.0f,1.0f,1.0f};
        glMaterialfv(GL_FRONT,GL_SPECULAR,matSpec);
        glMaterialf(GL_FRONT,GL_SHININESS,200);
        terreno.desenhar();
        agua.roda_agua(&terreno);

        if(desenhaAgua3d)
        agua.desenharDensidade();
    }else{

        glClear(GL_COLOR_BUFFER_BIT);

       agua.desenharDensidade2D();
       agua.roda_agua(&terreno);
    }

    erosao.roda_erosao();

}

void quit(){
    //saveData(&terreno);
    saveDataAndClose(&terreno);
    terreno.~Terreno();
    shallowWater.~ShallowWater();
    erosao.~Erosao();
}
int it = 0;
void display(){

    //codigoErosaoHidraulica1();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        if(panorama)
        {
            glTranslatef(0,0,-16.0);
            //glRotated(20.0,1.0,0.0,0.0);
            glRotated(60.0,1.0,0.0,0.0);
            glRotated(x,0.0,1.0,0.0);
            //x++;
        }else{
        //x++;
            glTranslatef(0,0,-10.0);
            glRotated(90.0,1.0,0.0,0.0);
        //glRotated(x,0.0,1.0,0.0);
        }

        if(rain)
        {
            shallowWater.fixedFont(3,3,50.0f);
            //shallowWater.rainDrop(500.0f);
        }

        //shallowWater.shallowWaterStep();
        //erosao.transportaSedimento();
        ///erosao.roda_erosao();
        //erosao.aeolianErosion();

        if(desenhaAgua3d)
        shallowWater.desenhar();
        terreno.desenhar();
        //erosao.desenha_sedimento();


    glFlush();
    glutSwapBuffers();
}

void pre_display_2D(){

    //glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);

    //glViewport(0,0,WIDTH,HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-10,10,-10,10);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void pre_display_3D(){


    glEnable(GL_DEPTH_TEST);

    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    float glGlobalAmb[] = {0.5f,0.5f,1.5f,0.5f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,glGlobalAmb);

    glEnable(GL_LIGHT0);
    float light0[4][4] = {

        {0.1f,0.1f,0.1f,1.0f}, // ambient
        {0.8f,0.8f,0.8f,1.0f}, // diffuse
        {1.0f,1.0f,1.0f,1.0f}, // specular
        {0.0f,0.0f,0.0f,1.0f} // position
    };
    glLightfv(GL_LIGHT0,GL_AMBIENT,&light0[0][0]);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,&light0[1][0]);
    glLightfv(GL_LIGHT0,GL_SPECULAR,&light0[2][0]);
    glLightfv(GL_LIGHT0,GL_POSITION,&light0[3][0]);

    //glViewport(0,0,WIDTH,HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60,1,1.0,50.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void reshape(int w, int h){
    glViewport(0,0,w,h);

    if(visual3D)
        pre_display_3D();
    else
        pre_display_2D();
}

void timer(int i){
    glutPostRedisplay();
    glutTimerFunc(1000/60,timer,0);
}

void keyboardFunction(unsigned char Key, int x, int y){

    switch(Key){
        case 'r':
            if(panorama)
            {
                //pre_display_2D();
                panorama = false;
            }else{
                //pre_display_3D();
                panorama = true;
            }
        break;
        case 'e':
            if(desenhaAgua3d)
            {
                desenhaAgua3d = false;
            }else{
                desenhaAgua3d = true;
            }
        break;
        case 'u':
            if(rain)
            {
                rain = false;
            }else{
                rain = true;
            }
        break;
        case 'q':
            agua.setDens(10,agua.N/2,agua.N/2);
        break;
        case 'd':
            agua.setU(10,agua.N/2,agua.N/2);
        break;
        case 'a':
            agua.setU(-10,agua.N/2,agua.N/2);
        break;
        case 'w':
            agua.setV(10,agua.N/2,agua.N/2);
        break;
        case 's':
            agua.setV(-10,agua.N/2,agua.N/2);
        break;
        case 'x':
            shallowWater.evaporateAllWater();
        break;
        case 'z':
            quit();
            exit(0);
        break;
    }

}
