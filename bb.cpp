#include "ibex.h"
#include "vibes.cpp"
#include <fstream>
#include <math.h> /*sqrt*/

#define __PREC__ 1e-11
#define __METH__ RK4
#define DELTA_T 0.1
#define EPS 1e-6
using namespace ibex;
using namespace std;

//Fonction d'aide pour afficher les résultats de la simulation
//TODO: La réécrire pour éviter l'erreur de segmentation
//void plot_simu(simulation sim)
//{
//    std::list<solution_g>::iterator iterator_list;
//    for(iterator_list=sim.list_solution_g.begin();iterator_list!=sim.list_solution_g.end();iterator_list++)
//    {
//      if (iterator_list->box_j1->diam().max() < 2)
//      {
//	vibes::drawBox(iterator_list->box_j1->operator[](0).lb(), iterator_list->box_j1->operator[](0).ub(), iterator_list->box_j1->operator[](1).lb(), iterator_list->box_j1->operator[](1).ub(), "k[g]");
//	std::cout << iterator_list->box_j1->operator[](0) << std::endl;
//      }
//    }
//
//}

//Créé un vecteur depuis l'intervalle [x] 
//constitué des scalaires x_min et x_max 
//À évaluer par la garde pour déterminer si le sol est traversé
IntervalVector interval_trans (Interval i){
  double _x[2][2]={{i.lb(),i.lb()},{i.ub(),i.ub()}};
  IntervalVector uplow (2,_x);
  return uplow;
}

//Imprime un intervalle
void print_interval(Interval i){
    printf("[%f,%f]",i.lb(),i.ub());
}

//Calcul du temps exact d'intersection avec le sol
//Entrée: conditions initiales en position et en vitesse verticales
double exact_intersect_time(double y0, double vy0){
    if (y0 < EPS){
        return (vy0/9.81);
    }
    else{
        return ( (-vy0 - sqrt(pow(vy0,2)+4*9.81*y0)) / (-2*9.81) );
    }
}

//Récupère les derniers résultats au dessus de la garde dans une simulation
//garantie: x, y, vx, vy, t
IntervalVector dernier_res_valide(simulation simu){
    Interval ytmp;
    list<solution_g>::iterator iterator_list;

    iterator_list=simu.list_solution_g.end();
    while(iterator_list!=simu.list_solution_g.begin()){
        ytmp = iterator_list->box_j1->operator[](1);
        if (ytmp.lb()<0){
            if (ytmp.ub()<0){
                iterator_list--; //Si on est en dessous du sol avec certitude, on passe au résultat précédent
            }
            else{
                IntervalVector res(5);
                res[0]=iterator_list->box_j1->operator[](0); //box containing x
                res[1]=iterator_list->box_j1->operator[](1); //box containing y
                res[2]=iterator_list->box_j1->operator[](2); //box containing vx
                res[3]=iterator_list->box_j1->operator[](3); //box containing vy
                res[4]=iterator_list->box_j1->operator[](4); //box containing t
                return res;
                break;
            }
        }
    }
}
int main() {
  
  double t_cross_exact = 0.0;
  double duration = 0.0;
  //**********
  //conditions initiales
  //***********
  //
  //Vitesse
  Interval Vx0(5.0,5.0);
  Interval Vy0(0.0,0.0);
  //Position
  Interval Px0(0.0,0.0);
  Interval Py0(100.0,100.0);
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
  

  
  //**********
  //dynamique
  //**********
  Variable y(5); 
  Function deriv = Function (y, Return(y[2],y[3],Interval(0.0,0.0),-g,Interval(1.0,1.0)) );

  //TODO:
  //**********
  //garde
  //***********
  Variable x1;
  Variable y1;
  Function guard = Function (x1,y1,x1*y1);
  //Nous donne un intervalle produit 
  Interval z = guard.eval(interval_trans(Py0));
  print_interval(z);

  //TODO:
  //**********
  //Résolution du problème
  //***********
  //Définition du problème
  ivp_ode problem = ivp_ode (deriv, 0.0, State_init);
  
  //Calcul du temps d'intersection exact avec le sol
  t_cross_exact = exact_intersect_time(State_init[1].lb(),State_init[3].lb());
  duration = t_cross_exact+DELTA_T;
  //Construction de la simulation, lancement
  simulation simu = simulation (&problem, duration, __METH__, __PREC__);
  simu.run_simulation();
  //Récupère les derniers résultats au dessus de la garde de la simulation courante
  //TODO
  //Bissection en temps
  //bissected_time = Bissect(t);
  //Définition d'un nouveau problème avec les conditions initiales bissectées
  //State_init[0] = Px0; // x
  //State_init[1] = Py0; // y
  //State_init[2] = Vx0; // vx
  //State_init[3] = Vy0; // vy
  //State_init[4] = T0;  // t
 
  // simu.export2d_yn("export-vitesse.txt", 2,3);
 // simu.export2d_yn("export-position.txt", 0,1);
 // simu.export1d_yn("export-time.txt",4);
  //TODO:
  //Dessin
  //vibes::beginDrawing ();
  //vibes::newFigure("plot");
  //plot_simu(simu);
  return 0;
}
