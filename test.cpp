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
//Function f(t,-9.81*pow(t,2) + 100); //Notre dynamique
Function ft(t, -9.81*pow(t,2)/2+100);
Function fy(y, y);
NumConstraint constraintOnT(ft,GEQ); //Contrainte: y(t) >= 0 
NumConstraint constraintOnY(fy,GEQ); //Contrainte: y(t) >= 0 

Variable ourVar(2); //(t,y)
Function fGlobal(ourVar, Return(-9.81*pow(ourVar[0],2)/2+100,ourVar[1]));
NumConstraint global(fGlobal,GEQ);

cout << "Intervalle test avant contraction: " << test << endl;
//Contraction
//CtcFwdBwd ctcOnT(constraintOnT);
//ctcOnT.contract(test);
//CtcFwdBwd ctcGlobal(global);
//ctcGlobal.contract(test0);
//cout << "Intervalle contracté T fwd-bwd" << test << endl;
// // Result : Intervalle contracté T fwd-bwd([4.065, 4.51524] ; [-20.914, 18.9489])
// cout << "Intervalle contracté global fwd-bwd" << test0 << endl;
// // Result : Intervalle contracté global fwd-bwd([4.065, 4.51524] ; [0, 18.9489])

CtcFwdBwd ctcGlobal(global);
CtcFixPoint fp(ctcGlobal, 1e-03);
fp.contract(test);
cout << "Intervalle contracté fwd-bwd et fixpoint" << test << endl;
// Result : Intervalle contracté fwd-bwd et fixpoint([4.065, 4.51524] ; [0, 18.9489])
// Remarque : aucun changement par rapport à seulement fwd-bwd --> quel est le meilleur contracteur à utiliser?? 

// // TODO: algo qui bissecte et contracte sur chaque intervalle ==> trouver des conditions sur x,y,vx,vy
// t : 0 <= -9.81/2 * T^2 + vY0 T + y0 <= epsilon
// x : min(vXo * T + xO ) <= X <= max(vX0 * T + x0 )
// y : 0 <= y <= epsilon
// vx : encadré par vx0 min et max?
// vy : encadré par min/max ( -9.81 t + vy0)
// REMARQUE : problème dépendance des contraintes => d'abord contracter sur t et y puis contracter sur les autres dimensions?

//Valeur initiales pour les fonctions :
//Vitesse
  Interval Vx0(5.0,5.0);
  Interval Vy0(0.0,0.0);
  //Position
  Interval Px0(0.0,0.0);
  Interval Py0(100.0,100.0);
  //Temps
  Interval T0(0.0,0.0);
  
  double Eps = 1e-6;
//STEP 0 : fonctions min et max qui permettent de trouver le min et le max sur des calculs d'intervalles? --> utilisation lb et ub


//STEP 1A : contracter sur T et Y (conditions >= 0 )
IntervalVector contractTY(2);
contractTY[0]=t_test;
contractTY[1]=y_test;
cout << "Intervalle TY avant contraction: " << contractTY << endl;
Variable ourVarTY(2); //(t,y)
double alpha = 0;
Function fA_TY(ourVarTY, Return(-9.81*pow(ourVarTY[0],2)/2+100 -alpha ,ourVarTY[1]- alpha));
NumConstraint A_TY(fA_TY,GEQ);
CtcFwdBwd A_ctcTY(A_TY);
CtcFixPoint A_fpTY(A_ctcTY, 1e-03);
A_fpTY.contract(contractTY);
cout << "Intervalle TY contracté (A) fwd-bwd et fixpoint" << contractTY<< endl;

// STEP 1B : contracter sur T et Y (conditions <= eps )
alpha = Eps;
Function fB_TY(ourVarTY, Return(-9.81*pow(ourVarTY[0],2)/2+100 -alpha ,ourVarTY[1]- alpha));
NumConstraint B_TY(fB_TY,LEQ);
CtcFwdBwd B_ctcTY(B_TY);
CtcFixPoint B_fpTY(B_ctcTY, 1e-03);
B_fpTY.contract(contractTY);
cout << "Intervalle TY contracté (B) fwd-bwd et fixpoint" << contractTY<< endl;


// // PROBLEME POUR CETTE CONTRACTION -> A MODIFIER
//STEP 2A : contracter X, Vx, Vy grâce à la contraction de TY (conditions >=0)
IntervalVector contractXVxVy(3);
contractXVxVy[0]=x_test;
contractXVxVy[1]=vx_test;
contractXVxVy[2]=vy_test;
cout << "Intervalle XVxVy avant contraction: " << contractXVxVy << endl;
Variable ourVarXVxVy(3); //(x,vx,vy)
Function A_fXVxVy(ourVarXVxVy, Return( ourVarXVxVy[0]-(Vx0*contractTY[0] + Px0).lb() , ourVarXVxVy[1]-Vx0.lb() , ourVarXVxVy[2]-(Vy0*contractTY[0]+Py0).lb() ) );
NumConstraint A_XVxVy(A_fXVxVy,GEQ);
CtcFwdBwd A_ctcXVxVy(A_XVxVy);
CtcFixPoint A_fpXVxVy(A_ctcXVxVy, 1e-03);
A_fpXVxVy.contract(contractXVxVy);
cout << "Intervalle XVxVy contracté (A) fwd-bwd et fixpoint" << contractXVxVy << endl;

// STEP 2B : contracter X, Vx, Vy grâce à la contraction de TY (conditions <=0)


//Version contraction des 5 en même temps => condition >=0
IntervalVector contractAll(5);
contractAll[0]=t_test;
contractAll[1]=y_test;
contractAll[2]=x_test;
contractAll[3]=vy_test;
contractAll[4]=vx_test;

cout << "Intervalle complet avant contraction : " << contractAll <<endl;
Variable all(5); //t,y,x,vy,vx
alpha = 0;
//Function A_f_const(all, Return( -9.81*pow(all[0],2)/2+100 -alpha , all[1]- alpha, all[2] - Vx0.lb()*all[0] - Px0.lb() , all[3] - Vy0.lb()*all[0] - Py0.lb() , all[4] - Vx0.lb() )); 
Function A_f_const(all, Return( -9.81*pow(all[0],2)/2+100 -alpha , all[1]- alpha, all[2] - Vx0.lb()*all[0] - Px0.lb() , all[3] - (-9.81*all[0] + Vy0.lb()) , all[4] - Vx0.lb() )); 
NumConstraint A_allConst(A_f_const,GEQ);
CtcFwdBwd A_ctc_All(A_allConst);
CtcFixPoint A_fpAll(A_ctc_All, 1e-5);
A_fpAll.contract(contractAll);
cout << "Intervalle complet contracté (A) " << contractAll << endl;


//Version contraction des 5 en même temps => condition <= eps
alpha = Eps;
//Function A_f_const(all, Return( -9.81*pow(all[0],2)/2+100 -alpha , all[1]- alpha, all[2] - Vx0.lb()*all[0] - Px0.lb() , all[3] - Vy0.lb()*all[0] - Py0.lb() , all[4] - Vx0.lb() )); 
Function B_f_const(all, Return( -9.81*pow(all[0],2)/2+100 -alpha , all[1]- alpha, all[2] - Vx0.ub()*all[0] - Px0.ub() , all[3] - (-9.81*all[0] + Vy0.ub()) , all[4] - Vx0.ub() )); 
NumConstraint B_allConst(B_f_const,LEQ);
CtcFwdBwd B_ctc_All(B_allConst);
CtcFixPoint B_fpAll(B_ctc_All, 1e-7);
B_fpAll.contract(contractAll);
cout << "Intervalle complet contracté (A) " << contractAll << endl;

//Besoin de bissecter les intervalles pour récupérer les valeurs de y > 0
//Question: est-ce qu'on ne peut pas brutalement retirer la partie y>0?

}
