#include "ibex.h"
#include "vibes.cpp"
#include <fstream>

#define __PREC__ 1e-6
#define __METH__ RK4
#define __DURATION__ 50.0

using namespace ibex;
using namespace std;

//Fonction d'aide pour afficher les résultats de la simulation

void plot_simu(simulation sim)
{
    std::list<solution_g>::iterator iterator_list;
    for(iterator_list=sim.list_solution_g.begin();iterator_list!=sim.list_solution_g.end();iterator_list++)
    {
      if (iterator_list->box_j1->diam().max() < 2)
      {
	vibes::drawBox(iterator_list->box_j1->operator[](0).lb(), iterator_list->box_j1->operator[](0).ub(), iterator_list->box_j1->operator[](1).lb(), iterator_list->box_j1->operator[](1).ub(), "k[g]");
	std::cout << iterator_list->box_j1->operator[](0) << std::endl;
      }
    }

}

int main() {
  
  vibes::beginDrawing ();
  vibes::newFigure("plot");

  //TODO: 
  
  //**********
  //conditions initiales
  //***********
  //IntervalVector yinit(2);
  //yinit[0] = Interval(0.0,0.0);
  //yinit[1] = Interval(5.0,5.0);
  //double g = -9.81; //gravité
  

  
  //**********
  //dynamique
  //**********
  //Variable y(2);
  //Function ydot = Function (g, y, Return (-g,
  //                                        y[1]) );


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
  //ivp_ode problem = ivp_ode (ydot, 0.0, yinit);

  //Construction de la simulation, lancement
  //simulation simu = simulation (&problem, __DURATION__, __METH__, __PREC__);
  //simu.run_simulation();

  //TODO:
  //Dessin
  plot_simu(simu);
  return 0;
}
