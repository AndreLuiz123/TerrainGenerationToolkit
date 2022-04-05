#ifndef TERRENO_H_INCLUDED
#define TERRENO_H_INCLUDED

#include <vector>

using namespace std;

typedef struct {

    float x;
    float y;
    float z;

}Vetor;

class Terreno{

private:
    int N=8;

public:
    vector<Vetor> alturas ;//= vector<Vetor>(N*N);
    vector<Vetor> alturas0 ;//= vector<Vetor>(N*N);
    vector<float> alturasY ;//= vector<Vetor>(N*N);
    vector<float> alturasY0 ;//= vector<Vetor>(N*N);

    Terreno();
    Terreno(int N, int numberHills=0);
    void desenhar();

    void inicializaAlturas();

    void hillAlgorithm(int numberHills);
    void faulting(int iterations);
    void valueNoise(vector<float> &values, unsigned seed=2001);
    void valueNoise2d(vector<float> &values, unsigned seed=2001);
    float coslerp(float a,float b,float alpha);
    float lerp(float a,float b,float alpha);
    float noise(vector<float> &values, float t);
    float noise2d(vector<float> &values, float tx, float ty);
    void perlinNoise();
    void diamond(int iterations);
    void diamond_squareStep(int chunk, float dHeight, float r);
    void diamond_diamondStep(int chunk, float dHeight, float r);
    void recursiveDiamond(int iterations, int x0,int y0, int h, float dHeight);
    void faultingDefaultFunction(int i, int j, float a, float b, float c, float displacement);
    void faultingSinFunction(int i, int j, float a, float b, float c, float displacement);
    void faultingCosFunction(int i, int j, float a, float b, float c, float displacement);
    void normalizeHeight(float maxH, float minH);
    float findMaxValue(bool maxOrMin);
    void salvarTerreno();
    void cubo();
    void zeraTerreno();
    void terrenoTeste();
    void variosCubos();
    void cubosUnidos();
    void terrenao();
    ~Terreno();

};

#endif // TERRENO_H_INCLUDED
