#include "ibex.h"
#include "vibes.cpp"
#include <fstream>
#include <math.h> /*sqrt*/
#include <algorithm> /*min max*/

#define __PREC__ 1e-11
#define __METH__ RK4
//#define DELTA_T 1.
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

// // INUTILE --> intégré dans le main pour question de seg fault
//Récupère les derniers résultats au dessus de la garde dans une simulation
//garantie: x, y, vx, vy, t
IntervalVector dernier_res_valide(simulation simu, list<solution_g>::iterator iterator_list){
//    Interval ytmp(0.,0.);
//    list<solution_g>::iterator mem, iterator_list;
    for(iterator_list=simu.list_solution_g.begin();iterator_list!=simu.list_solution_g.end();iterator_list++){
//        ytmp = iterator_list->box_j1->operator[](1);
//        print_interval(ytmp);
//        if (ytmp.lb()>0){
//            mem = iterator_list; //Si on est en dessous du sol avec certitude, on passe au résultat précédent
//        }
//        else{
//            if (ytmp.ub()>0){
//                //Si l'intervalle a un min négatif et un max positif, on
//                //renvoie le 
//                IntervalVector res(5);
//                res[0]=iterator_list->box_j1->operator[](0); //box containing x
//                res[1]=iterator_list->box_j1->operator[](1); //box containing y
//                res[2]=iterator_list->box_j1->operator[](2); //box containing vx
//                res[3]=iterator_list->box_j1->operator[](3); //box containing vy
//                res[4]=iterator_list->box_j1->operator[](4); //box containing t
//                return res;
//                break;
//            }
//            else{
//            //Si on a une boîte au dessus de 0 et une boîte en dessous de 0, on
//            //fait une union brutale entre les intervalles
//            print_interval(iterator_list->box_j1->operator[](0));
//            print_interval(mem->box_j1->operator[](0));
//            IntervalVector res(5);
//            
//            res[0]=union_interval(iterator_list->box_j1->operator[](0),mem->box_j1->operator[](0)); //box containing x from now
//            print_interval(res[0]);
//            res[1]=union_interval(iterator_list->box_j1->operator[](1),mem->box_j1->operator[](1)); //box containing y from now
//            res[2]=union_interval(iterator_list->box_j1->operator[](2),mem->box_j1->operator[](2)); //box containing vx from now
//            res[3]=union_interval(iterator_list->box_j1->operator[](3),mem->box_j1->operator[](3)); //box containing vy from now
//            res[4]=union_interval(iterator_list->box_j1->operator[](4),mem->box_j1->operator[](4)); //box containing t from now
//            return res;
//            break;
//            }
//        }
// //FIN INUTILE
    }
}
int main() {
  
  double t_cross_exact = 0.0;
  double duration = 0.0;
  double Eps = 1e-3;
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
  //Intervalle à bissecter
  IntervalVector ToBB(5);
  ToBB[0]=Px0;
  ToBB[1]=Px0;
  ToBB[2]=Px0;
  ToBB[3]=Px0;
  ToBB[4]=Px0;
  
  //**********
  //dynamique
  //**********
  Variable y(5); 
  Function deriv = Function (y, Return(y[2],y[3],Interval(0.0,0.0),-g,Interval(1.0,1.0)) );

  ////TODO:for later maybe
  ////**********
  ////garde
  ////***********
  //Variable x1;
  //Variable y1;
  //Function guard = Function (x1,y1,x1*y1);
  ////Nous donne un intervalle produit 
  //Interval z = guard.eval(interval_trans(Py0));
  //print_interval(z);

  //**********
  //Résolution du problème
  //***********
  //Définition du problème initial
  ivp_ode problem = ivp_ode (deriv, 0.0, State_init);
  
  //Calcul du temps d'intersection exact avec le sol
  t_cross_exact = exact_intersect_time(State_init[1].mid(),State_init[3].mid());
  double delta_T = exact_intersect_time(State_init[1].diam(),State_init[3].diam());
  duration = t_cross_exact+delta_T;
  cout << "T simulation " << delta_T << endl;
  cout << "T exact d'intersection " << t_cross_exact <<endl;
  cout << "Y exact d'intersection " << -9.81 * pow(t_cross_exact,2)/2 + Vy0.mid() * t_cross_exact + Py0.mid() << endl;
  cout << "X exact d'intersection " << Vx0.mid() * t_cross_exact + Px0.mid() << endl;
  cout << "Vy exact d'intersection " << -9.81 *t_cross_exact + Vy0.mid()<<endl;
  cout << "Vx exact d'intersection " << Vx0.mid() <<endl;
  //Construction de la simulation, lancement
  simulation simu = simulation (&problem, duration, __METH__, __PREC__);
  simu.run_simulation();
  simu.export2d_yn("export-vitesse.txt", 2,3);
  simu.export2d_yn("export-position.txt", 0,1);
  simu.export2d_yn("export-time.txt",0,4);
  //Récupère les derniers résultats au dessus de la garde de la simulation courante
  list<solution_g>::iterator it,mem;
  Interval ytmp;
  for(it=simu.list_solution_g.begin();it!=simu.list_solution_g.end();it++){
    ytmp = it->box_j1->operator[](1);
    if (ytmp.lb()>0){
        mem = it; //Si on est eu dessus du sol avec certitude, on passe au résultat suivant
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
            print_interval(ToBB[0]);
            ToBB[1]=union_interval(it->box_j1->operator[](1),mem->box_j1->operator[](1));
            ToBB[2]=union_interval(it->box_j1->operator[](2),mem->box_j1->operator[](2));
            ToBB[3]=union_interval(it->box_j1->operator[](3),mem->box_j1->operator[](3));
            ToBB[4]=union_interval(it->box_j1->operator[](4),mem->box_j1->operator[](4));
            break;
        }
    }
  }
  printf("Intervalle ToBB : intersection \n ");
  for (int i=0;i<5;i++){
    printf("\n");
    print_interval(ToBB[i]);
  }
  
  
  
  //Version contraction des 5 en même temps => condition y>=0 && vx >= min(vx(t)) &&  x>= min(x(t)) && vy >= min(vy(t))
  cout << "Intervalle complet avant contraction : " << ToBB <<endl;
  Variable all(5); //x,y,vx,vy,t
  double alpha = 0;
  Function A_f_const(all, Return( all[0] - Vx0.lb()*all[4] - Px0.lb() ,  all[1]- alpha,  all[2] - Vx0.lb(), all[3] - (-9.81*all[4] + Vy0.lb()) , -9.81*pow(all[4],2)/2+100 -alpha )); 
  NumConstraint A_allConst(A_f_const,GEQ);
  CtcFwdBwd A_ctc_All(A_allConst);
  CtcFixPoint A_fpAll(A_ctc_All, 1e-7);
  A_fpAll.contract(ToBB);
  cout << "Intervalle complet contracté (contrainte y>=0) " << ToBB << endl;


  //Version contraction des 5 en même temps => condition y<=eps && vx <= max(vx(t)) &&  x<= max(x(t)) && vy <= max(vy(t))
  alpha = Eps;
  Function B_f_const(all, Return( all[0] - Vx0.ub()*all[4] - Px0.ub() ,  all[1]- alpha,  all[2] - Vx0.ub(), all[3] - (-9.81*all[4] + Vy0.ub()) , -9.81*pow(all[4],2)/2+100 -alpha )); 
  NumConstraint B_allConst(B_f_const,LEQ);
  CtcFwdBwd B_ctc_All(B_allConst);
  CtcFixPoint B_fpAll(B_ctc_All, 1e-4);
  B_fpAll.contract(ToBB);
  cout << "Intervalle complet contracté (contrainte y<=epsilon) " << ToBB << endl;

  
  
  
  //TODO: Tester des stratégies de contraction sur les différents intervalles
  //Définition d'un nouveau problème avec les conditions initiales bissectées
  //State_init[0] = Px0; // x
  //State_init[1] = Py0; // y
  //State_init[2] = Vx0; // vx
  //State_init[3] = Vy0; // vy
  //State_init[4] = T0;  // t
  //problem = problem(...,State_init);
  //Relancer la simulation
  //
 
  //TODO:
  //Dessin
  //vibes::beginDrawing ();
  //vibes::newFigure("plot");
  //plot_simu(simu);
  return 0;
}
