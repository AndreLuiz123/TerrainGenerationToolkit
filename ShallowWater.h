#ifndef SHALLOWWATER_H_INCLUDED
#define SHALLOWWATER_H_INCLUDED

#include "Terreno.h"

typedef struct {

    float L=0;
    float R=0;
    float T=0;
    float B=0;
    float K=0;

}Flow;

class ShallowWater{

private:

public:
    int N;
    vector<Vetor> alturaAgua;
    vector<Vetor> alturaAguaVelha;
    vector<Flow> flow;
    Terreno *terreno;
    float dt=185.1f;
    float Ix, Iy;
    vector<float> u, v;
    float kr = 1;

    ShallowWater(int N, Terreno *terreno, float dt);

    void inicializaAlturas();
    int inicializaVelocidadeAgua();
    void free_data();
    void atualizaAlturas(int i, int j);
    void desenhar();
    void waterIncrement(int i, int j, float r);
    void waterEvaporation();
    void rainDrop(float r=0);
    void fixedFont(int i, int j, float r=0);
    void riverSource();
    void inflowFlux(int i, int j);
    void outflowFlux(int i, int j);
    void newVolume(int i, int j);
    void shallowWaterStep();
    void atualizaAlturaVelha();
    void evaporateAllWater();
    void trataBordas();
    ~ShallowWater();

};

#endif // SHALLOWWATER_H_INCLUDED
