#include "ibex.h"
#include "vibes.cpp"
#include <fstream>

#define __PREC__ 1e-11
#define __METH__ RK4
#define __DURATION__ 10.0

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

int main() {
  


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
  //Vecteur d'état initial
  IntervalVector State_init(4);
  State_init[0] = Px0; // x
  State_init[1] = Py0; // y
  State_init[2] = Vx0; // vx
  State_init[3] = Vy0; // vy
  //Constante
  Interval g(9.81,9.81); //gravité
  

  
  //**********
  //dynamique
  //**********
  Variable y(4); 
  Function deriv = Function (y, Return(y[2],y[3],Interval(0.0,0.0),-g) );

  //TODO:
  //**********
  //garde
  //***********
  //Function g = Function (...,Return(...))
  

  //TODO:
  //**********
  //définition du problème
  //***********
  //Définition du problème
  ivp_ode problem = ivp_ode (deriv, 0.0, State_init);

  //Construction de la simulation, lancement
  simulation simu = simulation (&problem, __DURATION__, __METH__, __PREC__);
  simu.run_simulation();
  simu.export2d_yn("export-vitesse.txt", 2,3);
  simu.export2d_yn("export-position.txt", 0,1);
  //TODO:
  //Dessin
  //vibes::beginDrawing ();
  //vibes::newFigure("plot");
  //plot_simu(simu);
  return 0;
}
