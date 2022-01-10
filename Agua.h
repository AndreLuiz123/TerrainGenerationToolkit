#ifndef AGUA_H_INCLUDED
#define AGUA_H_INCLUDED

#include <vector>
#include "Terreno.h"
using namespace std;

class Agua{

private:
public:
    int N=8;
    float * u, * v, * u_prev, * v_prev;
    float * dens, * dens_prev;
    float dt, diff=0.5f, visc;
    float force, source;
    int dvel;
    int dist_terreno;
    Agua();
    void desenharDensidade();
    void desenharDensidade2D();
    void desenharVelocidade();
    void add_source(int N, float * x, float * s, float dt);
    void advect(int N, int b, float * d, float * d0, float * u, float * v, float dt);
    void project(int N, float * u, float * v, float * p, float * div);
    void diffuse(int N, int b, float * x, float * x0, float diff, float dt);
    void lin_solve(int N, int b, float * x, float * x0, float a, float c);
    void set_bnd(int N, int b, float * x);

    void advectDens(int N, int b, vector<Vetor> &d, vector<Vetor> &d0, float * u, float * v, float dt);
    void diffuseDens(int N, int b, vector<Vetor> &x, vector<Vetor> &x0, float diff, float dt);
    void lin_solve_Dens(int N, int b, vector<Vetor> &x, vector<Vetor> &x0, float a, float c);
    void set_bnd_Dens(int N, int b, vector<Vetor> &x);

    void vel_step(int N, float * u, float * v, float * u0, float * v0, float visc, float dt);
    void dens_step(int N, vector<Vetor> &x, vector<Vetor> &x0, float * u, float * v, float diff, float dt);
    void roda_agua(Terreno *terreno);

    void setDens(float value, int i, int j);
    void setU(float value, int i, int j);
    void setV(float value, int i, int j);

    void free_data();
    void clear_data();
    int allocate_data();
    ~Agua();
};


#endif // AGUA_H_INCLUDED
