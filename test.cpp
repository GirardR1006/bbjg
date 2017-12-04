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
Interval vy_test(-48.707,-39.877);
Interval vx_test(5.0,5.0);
IntervalVector test(2);
test[0]=t_test;
test[1]=y_test;
IntervalVector test0(2);
test0=test;

Variable y;
Variable x;
Variable t;

Function ft(t, -9.81*pow(t,2)/2+100); //Dynamique
NumConstraint constraintOnT(ft,GEQ); //Contrainte: y(t) >= 0 

// Contraintes 2 dimensions
Variable ourVar(2); //(t,y)
Function fGlobal(ourVar, Return(-9.81*pow(ourVar[0],2)/2+100,ourVar[1])); //retourne [ f(t) ; y ]
NumConstraint global(fGlobal,GEQ);

cout << "Intervalle test avant contraction: " << test << endl;
// // Contraction 1D
//CtcFwdBwd ctcOnT(constraintOnT);
//ctcOnT.contract(test);
  
// // Contraction 2D
//CtcFwdBwd ctcGlobal(global);
//ctcGlobal.contract(test0);
//cout << "Intervalle contracté T fwd-bwd" << test << endl;
// // Result : Intervalle contracté T fwd-bwd([4.065, 4.51524] ; [-20.914, 18.9489])
// cout << "Intervalle contracté global fwd-bwd" << test0 << endl;
// // Result : Intervalle contracté global fwd-bwd([4.065, 4.51524] ; [0, 18.9489])

// // Contraction 2D composée
CtcFwdBwd ctcGlobal(global);
CtcFixPoint fp(ctcGlobal, 1e-03);
fp.contract(test);
cout << "Intervalle contracté fwd-bwd et fixpoint" << test << endl;
// Result : Intervalle contracté fwd-bwd et fixpoint([4.065, 4.51524] ; [0, 18.9489])
// Remarque : aucun changement par rapport à seulement fwd-bwd --> quel est le meilleur contracteur à utiliser?? 

//Besoin de bissecter les intervalles pour récupérer les valeurs de y > 0
//Question: est-ce qu'on ne peut pas brutalement retirer la partie y>0?

}
