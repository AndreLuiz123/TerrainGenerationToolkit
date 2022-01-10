#ifndef EROSAO_H_INCLUDED
#define EROSAO_H_INCLUDED

#include "Terreno.h"
#include "ShallowWater.h"
#include "Agua.h"

class Erosao{

private:

public:
    int N;
    float kc;
    float kd;
    float ks;
    float kt;
    float ka;
    //Terreno *terreno;
    ShallowWater *agua;
    Terreno *terreno;
    Agua *vento;
    //Agua *agua2;
    vector<float> sedimento_suspenso;
    vector<float> sedimento_suspenso0;
    float dt;

    Erosao(int N, Terreno *terreno, ShallowWater *agua, float kc, float kd, float ks);

    void inicializaSedimentoSuspenso();
    void roda_erosao();
    void roda_erosao_termica();
    void deposita_dissolve();
    void dissolveSedimento(int i,int j,float ds);
    void depositaSedimento(int i,int j,float ds);
    void transportaSedimento();
    void transportaSedimento2();
    void guardaGradeSedimentoSuspenso();
    void desenha_sedimento();
    void thermalErosion();
    void aeolianErosion();

    ~Erosao();

};

#endif // EROSAO_H_INCLUDED
