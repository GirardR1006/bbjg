#include "ibex.h"
#include "vibes.cpp"
#include <fstream>
#include <math.h> //pow

#define __PREC__ 1e-11
#define __METH__ RK4
#define DELTA_T 1.
#define EPS 1e-6
using namespace ibex;
using namespace std;

int main(){
//Valeurs renvoyées par l'exécutable principal
Interval y_test(-20.914019,18.948873);
Interval x_test(20.325,24.925);
Interval t_test(4.065,4.965);
IntervalVector test(2);
test[0]=t_test;
test[1]=y_test;


Variable y;
Variable x;
Variable t;
Function f(t,-9.81*pow(t,2) + 100); //Notre dynamique
NumConstraint c(f,GEQ); //Contrainte: y(t) >= 0 


cout << "Intervalle test avant contraction: " << test << endl;
//Contraction
CtcFwdBwd ctc(c);
CtcFixPoint fp(ctc,1e-03);
fp.contract(test);
cout << "Intervalle test après contraction: " << test << endl;
//Pour l'instant, pas de résultat probant (j'ai un empty vector en sortie) 

//Besoin de bissecter les intervalles pour récupérer les valeurs de y > 0
//Question: est-ce qu'on ne peut pas brutalement retirer la partie y>0?
}
