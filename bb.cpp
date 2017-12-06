#include "ibex.h"
#include "vibes.cpp"
#include <fstream>
#include <math.h> /*sqrt*/
#include <algorithm> /*min max*/

#define __PREC__ 1e-11
#define __METH__ RK4
#define EPS 1e-6
using namespace ibex;
using namespace std;

//Fonction d'aide pour afficher les résultats de la simulation
//Les boîtes dessinées sont les valeurs y(t)
void plot_all(list<IntervalVector> l)
{
    for(list<IntervalVector>::iterator it = l.begin();it!=l.end();it++){
	vibes::drawBox(it->operator[](4).lb(),it->operator[](4).ub(),it->operator[](1).lb(),it->operator[](1).ub(), "[red]");
       // cout << *it << endl;	
    }

}
//Calcul du temps exact d'intersection avec le sol
//Entrée: conditions initiales en position et en vitesse verticales
double exact_intersect_time(double y0, double vy0){
    if (y0 < EPS){
        return (2*vy0/9.81);
    }
    else{
        return ( (vy0 + sqrt(pow(vy0,2)+2*9.81*y0)) / (9.81) );
    }
}

Interval union_interval (Interval i1, Interval i2){
    Interval u(min(i1.lb(),i2.lb()),max(i1.ub(),i2.ub()));
    return u;
}
//Renvoie les intervalles y0 qui vérifient les conditions de rentrée 
//dans le panier. 
//Paramètre d'entrée: int n: nombre de rebond
//Paramètre d'entrée: IntervalVector init: vecteur de condition initiales
//Paramètre d'entrée: Function f: dynamique
//list<Interval> simulate_launch(int n, IntervalVector init, Function f){
//    IntervalVector w_in = init; //Conditions initiales amenées à être modifiées 
//    IntervalVector ToBB = init; //Conditions initiales amenées à être bissectée
//    double alpha = 0;
//    double t_cross_exact = 0;
//    list<Interval> valid_y;
//    for(int i=0;i<=n;i++){
//        //Définition des contracteurs
//        //condition y>=0 && vx >= min(vx(t)) &&  x>= min(x(t)) && vy >= min(vy(t))
//        Variable all(5); //x,y,vx,vy,t
//        Function A_f_const(all, Return( all[0] - w_in[2].lb()*all[4] - w_in[0].lb() 
//                                      , all[1]- alpha
//                                      , all[2] - w_in[2].lb()
//                                      , all[3] - (-9.81*all[4] + w_in[3].lb()) 
//                                      , -9.81*pow(all[4],2)/2+w_in[3].ub()*all[4] + w_in[1].ub() -alpha ));  
//        NumConstraint A_allConst(A_f_const,GEQ);
//        CtcFwdBwd A_ctc_All(A_allConst);
//        CtcFixPoint A_fpAll(A_ctc_All, 1e-7);
//
//         //condition y<=eps && vx <= max(vx(t)) &&  x<= max(x(t)) && vy <= max(vy(t))
//        Function B_f_const(all, Return( all[0] - w_in[2].ub()*all[4] - w_in[0].ub() 
//                                      , all[1]- alpha
//                                      , all[2] - w_in[2].ub() 
//                                      , all[3] - (-9.81*all[4] + w_in[3].ub())  
//                                      , -9.81*pow(all[4],2)/2+Vy0.lb()*all[4] + w_in[1].lb() -alpha )); 
//        NumConstraint B_allConst(B_f_const,LEQ);
//        CtcFwdBwd B_ctc_All(B_allConst);
//        CtcFixPoint B_fpAll(B_ctc_All, 1e-4);
//        //Définition du problème
//        ivp_ode problem = ivp_ode (f, w_in[4].mid(), w_in);
//        //Calcul du temps d'intersection exact avec le sol
//        t_cross_exact = exact_intersect_time(w_in[1].mid(),w_in[3].mid());
//        double delta_T = exact_intersect_time(w_in[1].diam(),w_in[3].diam());
//        duration = t_cross_exact+delta_T;
//        simulation simu = simulation (&problem, duration, __METH__, __PREC__);
//        simu.run_simulation();
//        //Récupère les derniers résultats au dessus de la garde de la simulation courante
//        list<solution_g>::iterator it,mem;
//        Interval ytmp;
//        for(it=simu.list_solution_g.begin();it!=simu.list_solution_g.end();it++){
//            ytmp = it->box_j1->operator[](1);
//            if (ytmp.lb()>0){
//                mem = it; 
//            //Si on est eu dessus du sol avec certitude, on passe au résultat suivant
//            }
//            else{
//                if (ytmp.ub()>0){
//                    ToBB[0]=it->box_j1->operator[](0); // x
//                    ToBB[1]=it->box_j1->operator[](1); // y
//                    ToBB[2]=it->box_j1->operator[](2); // vx
//                    ToBB[3]=it->box_j1->operator[](3); // vy
//                    ToBB[4]=it->box_j1->operator[](4); // t
//                    break;
//                }
//                else{
//            //Si on a une boîte au dessus de 0 et une boîte en dessous de 0, 
//            //on fait une union brutale entre les intervalles
//                    ToBB[0]=union_interval(it->box_j1->operator[](0),mem->box_j1->operator[](0));
//                    ToBB[1]=union_interval(it->box_j1->operator[](1),mem->box_j1->operator[](1));
//                    ToBB[2]=union_interval(it->box_j1->operator[](2),mem->box_j1->operator[](2));
//                    ToBB[3]=union_interval(it->box_j1->operator[](3),mem->box_j1->operator[](3));
//                    ToBB[4]=union_interval(it->box_j1->operator[](4),mem->box_j1->operator[](4));
//                    break;
//                }
//             }
//         }
//   
//    A_fpAll.contract(ToBB);
//    alpha = Eps;
//    B_fpAll.contract(ToBB);
//    //Vérifier les contraintes sur x ici.
//    //if contraintes_x not verified then 
//        //diviser w_in en deux w_in1 et w_in2
//        //simulate_launch(n, w_in1, f)
//        //simulate_launch(n, w_in2, f)
//
//    w_int = ToBB;  //Nouvelles conditions initiales 
//    valid_y.push_back(ToBB[1]);
//    }
//    //Vérifier les contraintes sur panier ici.
//    //if contraintes_y not verified then 
//        //diviser w_in en deux w_in1 et w_in2
//        //simulate_launch(n, w_in1, f)
//        //simulate_launch(n, w_in2, f)
//}

int main() {
  
  double t_cross_exact = 0.0;
  double duration = 0.0;
  double Eps = 1e-2;
  //**********
  //conditions initiales
  //**********
  //Vitesse
  Interval Vx0(5.0,5.0);
  Interval Vy0(0.0,0.0);
  //Position
  Interval Px0(0.0,0.0);
  Interval Py0(1.6,2.8);
  //Temps
  Interval T0(0.0,0.0);
  //Vecteur d'état initial
  IntervalVector State_init(5);
  State_init[0] = Px0; // x
  State_init[1] = Py0; // y
  State_init[2] = Vx0; // vx
  State_init[3] = Vy0; // vy
  State_init[4] = T0;  // t
  //Constante
  Interval g(9.81,9.81); //gravité
  double c = 0.8;   //elasticité du ballon
  //Intervalle à bissecter
  IntervalVector ToBB(5);
  ToBB[0]=Px0; //Initialisé à Px0 mais ça pourrait être n'importe quoi, à changer
  ToBB[1]=Px0;
  ToBB[2]=Px0;
  ToBB[3]=Px0;
  ToBB[4]=Px0;
  
  //**********
  //dynamique
  //**********
  Variable y(5); 
  Function deriv = Function (y, Return(y[2],y[3],Interval(0.0,0.0),-g,Interval(1.0,1.0)) );
  
  //*********
  //dessin
  //*********
  list<IntervalVector> toPlot;
  

  //*********
  //Résolution du problème
  //*********
  //Définition du problème initial
  ivp_ode problem = ivp_ode (deriv, 0.0, State_init);

 

  //Calcul du temps d'intersection exact avec le sol
  t_cross_exact = exact_intersect_time(State_init[1].mid(),State_init[3].mid());
  double delta_T = exact_intersect_time(State_init[1].diam(),State_init[3].diam());
  duration = t_cross_exact+delta_T;
  //cout << "T simulation " << delta_T << endl;
  //cout << "T exact d'intersection " << t_cross_exact <<endl;
  //cout << "Y exact d'intersection " << -9.81 * pow(t_cross_exact,2)/2 + Vy0.mid() * t_cross_exact + Py0.mid() << endl;
  //cout << "X exact d'intersection " << Vx0.mid() * t_cross_exact + Px0.mid() << endl;
  //cout << "Vy exact d'intersection " << -9.81 *t_cross_exact + Vy0.mid()<<endl;
  //cout << "Vx exact d'intersection " << Vx0.mid() <<endl;
  //Construction de la simulation, lancement
  simulation simu = simulation (&problem, duration, __METH__, __PREC__);
  simu.run_simulation();
  //simu.export2d_yn("export-vitesse.txt", 2,3);
  //simu.export2d_yn("export-position.txt", 0,1);
  //simu.export2d_yn("export-time.txt",0,4);
  //Récupère les derniers résultats au dessus de la garde de la simulation courante
  list<solution_g>::iterator it,mem;
  Interval ytmp;
  for(it=simu.list_solution_g.begin();it!=simu.list_solution_g.end();it++){
    ytmp = it->box_j1->operator[](1);
    if (ytmp.lb()>0){
        mem = it; //Si on est eu dessus du sol avec certitude, on passe au résultat suivant
        toPlot.push_back(*it->box_j1);
    }
    else{
        if (ytmp.ub()>0){
            ToBB[0]=it->box_j1->operator[](0); // x
            ToBB[1]=it->box_j1->operator[](1); // y
            ToBB[2]=it->box_j1->operator[](2); // vx
            ToBB[3]=it->box_j1->operator[](3); // vy
            ToBB[4]=it->box_j1->operator[](4); // t
            break;
        }
        else{
        //Si on a une boîte au dessus de 0 et une boîte en dessous de 0, on
        //fait une union brutale entre les intervalles
            ToBB[0]=union_interval(it->box_j1->operator[](0),mem->box_j1->operator[](0));
            ToBB[1]=union_interval(it->box_j1->operator[](1),mem->box_j1->operator[](1));
            ToBB[2]=union_interval(it->box_j1->operator[](2),mem->box_j1->operator[](2));
            ToBB[3]=union_interval(it->box_j1->operator[](3),mem->box_j1->operator[](3));
            ToBB[4]=union_interval(it->box_j1->operator[](4),mem->box_j1->operator[](4));
            break;
        }
    }
  }
   
  //Version contraction des 5 en même temps => condition y>=0 && vx >= min(vx(t)) &&  x>= min(x(t)) && vy >= min(vy(t))
 // cout << "Intervalle complet avant contraction : " << ToBB <<endl;
  Variable all(5); //x,y,vx,vy,t
  double alpha = 0;
  Function A_f_const(all, Return( all[0] - Vx0.lb()*all[4] - Px0.lb() ,  all[1]- alpha,  all[2] - Vx0.lb(), all[3] - (-9.81*all[4] + Vy0.lb()) , -9.81*pow(all[4],2)/2+Vy0.ub()*all[4] + Py0.ub() -alpha )); 
  
  NumConstraint A_allConst(A_f_const,GEQ);
  CtcFwdBwd A_ctc_All(A_allConst);
  CtcFixPoint A_fpAll(A_ctc_All, 1e-7);
  A_fpAll.contract(ToBB);
  toPlot.push_back(ToBB);
 // cout << "Intervalle complet contracté (contrainte y>=0) " << ToBB << endl;

  //Version contraction des 5 en même temps => condition y<=eps && vx <= max(vx(t)) &&  x<= max(x(t)) && vy <= max(vy(t))
  alpha = Eps;
  Function B_f_const(all, Return(all[0] - Vx0.ub()*all[4] - Px0.ub() , all[1]- alpha, all[2] - Vx0.ub() ,all[3] - (-9.81*all[4] + Vy0.ub())  , -9.81*pow(all[4],2)/2+Vy0.lb()*all[4] + Py0.lb() -alpha )); 
  NumConstraint B_allConst(B_f_const,LEQ);
  CtcFwdBwd B_ctc_All(B_allConst);
  CtcFixPoint B_fpAll(B_ctc_All, 1e-4);
  B_fpAll.contract(ToBB);
 // cout << "Intervalle complet contracté (contrainte y<=epsilon) " << ToBB << endl;
  
  toPlot.push_back(ToBB);
  cout << "Voici ToBB: " << endl;
  cout << ToBB << endl;
  //Définition d'un nouveau problème avec les conditions initiales bissectées et
  //la nouvelle dynamique
  State_init[0] = ToBB[0];  // x
  State_init[1] = ToBB[1];  // y
  State_init[2] = c*ToBB[2];  // vx
  State_init[3] = (-c)*ToBB[3];  // vy, condition initiale négative pour simuler le rebond
  State_init[4] = ToBB[4];  // t
  ivp_ode nproblem = ivp_ode (deriv, ToBB[4].mid(), State_init);
  
  t_cross_exact = t_cross_exact + exact_intersect_time(State_init[1].mid(),State_init[3].mid());
  delta_T = exact_intersect_time(State_init[1].diam(),State_init[3].diam());
  duration = t_cross_exact+delta_T;
  simulation nsimu = simulation (&nproblem,t_cross_exact, __METH__, __PREC__);
  //Relancer la nouvelle simulation
  nsimu.run_simulation();
 
  for(it=nsimu.list_solution_g.begin();it!=nsimu.list_solution_g.end();it++){
      toPlot.push_back(*it->box_j1);
  }
  //TODO:
  //Dessin
  vibes::beginDrawing ();
  vibes::newFigure("Basket");
  //for(list<IntervalVector>::iterator it1 = toPlot.begin();it1!=toPlot.end();it1++){
  //      cout << "Vecteur:"  << endl;	
  //      cout << *it1 << endl;	
  //  }
  plot_all(toPlot);

  return 0;
}
