//TODO:
//*Le delta_t calculé ne permet parfois pas d'atteindre la solution (on arrive
//pas à générer un nouveau toBB)
//*Des conditions initiales négatives sur la vitesse verticale (équivalent du
//slam dunk) provoquent un segfault sur la seconde contraction
//*Dans le cas d'un x trop large, les contracteurs le transforment en empty
//vector, ce qui empêche au pas suivant de la simulation de se dérouler
//Cela apparaît pour des conditions initiales non-nulles sur Vy
//*Une erreur "bad::alloc quand on a deux assignations de toBB d'affilée.


#include "ibex.h"
#include "vibes.cpp"
#include <fstream>
#include <math.h> /*sqrt*/
#include <algorithm> /*min max*/
#include <tuple>

#define __PREC__ 1e-11
#define __METH__ RK4
#define EPS 1e-6
using namespace ibex;
using namespace std;

//Fonction d'aide pour afficher les résultats de la simulation
//Les boîtes dessinées sont les valeurs y(t)
void plot_all(list<IntervalVector> l){
    for(list<IntervalVector>::iterator it = l.begin();it!=l.end();it++){
	vibes::drawBox(it->operator[](4).lb(),it->operator[](4).ub(),it->operator[](1).lb(),it->operator[](1).ub(), "[red]");
//        cout <<"Intervalle dessiné: " << *it << endl;	
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
//Union d'intervalle
Interval union_interval (Interval i1, Interval i2){
    Interval u(min(i1.lb(),i2.lb()),max(i1.ub(),i2.ub()));
    return u;
}
//Division d'intervalle en deux
tuple<Interval,Interval> divide_interval(Interval i){
    Interval i1 = Interval(i.lb(),i.mid());
    Interval i2 = Interval(i.mid(),i.ub());
    tuple<Interval,Interval> two_divided_interval(i1,i2);
    return two_divided_interval;
}
//Fonction vérifiant si xballe.inf > 6.95
//Paramètre d'entrée: IntervalVector test: vecteur à tester
//Paramètre d'entrée: IntervalVector safe: vecteur de zone panier
bool is_behind_basket (IntervalVector test,IntervalVector safe){
    return (test[0].lb()>safe[0].ub());
}
//Fonction vérifiant si une simulation franchit bien la boîte délimitant le panier
//Paramètre d'entrée: IntervalVector test: vecteur à tester
//Paramètre d'entrée: IntervalVector safe: vecteur de zone panier
//On s'attend à ce que la dimension en vy du vecteur safe soit négative
bool has_crossed_basket (simulation* simu,IntervalVector objective){
    return(simu->has_crossed(objective));
}
//L'objectif ici est de récupérer les boîtes résultats de la simulation.
//Comme la simulation dure le temps nécessaire à la balle pour retomber,
//la dernière boîte (après bissection) représentera notre nouvelle condition
//initiale
//Paramètre d'entrée: simulation* simu: résultats de la simulation 
//Paramètre d'entrée: list<IntervalVector>* toPlot: boîtes à tracer
//Paramètre d'entrée: IntervalVector toBB: boîte à bissecter
IntervalVector get_toBB(simulation* simu, list<IntervalVector>* toPlot, IntervalVector toBB){ 
    Interval ytmp, vytmp;
    list<solution_g>::iterator it,mem; //itérateurs sur la liste des solutions et mémoire de la solution précédente
    mem=simu->list_solution_g.begin(); 
    for(it=simu->list_solution_g.begin();it!=simu->list_solution_g.end();it++){
        ytmp  = it->box_j1->operator[](1);
        vytmp = it->box_j1->operator[](3);
        cout << "Solution: " << *it->box_j1 << endl;
        //Si on est eu dessus du sol avec certitude, on passe à la prochaine itération
        //en sauvegardant la solution courante pour comparaison ultérieure
        if (ytmp.lb()>=0){
            mem = it; 
            toPlot->push_back(*it->box_j1);
        }
        else{
        //Si l'intervalle auquel appartient y contient 0 et que y décroît,
        //alors la balle rebondira dans l'intervalle de temps décrit par la boîte
            if (ytmp.ub()>0 && vytmp.lb()<=0){
                toBB=*it->box_j1; 
//                toPlot.push_back(*it->box_j1);
                cout << "New toBB is here! " << toBB << endl;
                return(toBB); 
            }
            else{
       //Si on a une boîte aux valeurs de y strictement négatives et que la
       //boîte précédente avait des valeurs de y strictements positives, alors
       //on fait l'union de ces deux intervalles. On vérifie que la vitesse
       //est décroissante pour ne pas prendre les boîtes ascendantes 
                if(vytmp.ub()<=0){
                    toBB[0]=union_interval(it->box_j1->operator[](0),mem->box_j1->operator[](0));
                    toBB[1]=union_interval(it->box_j1->operator[](1),mem->box_j1->operator[](1));
                    toBB[2]=union_interval(it->box_j1->operator[](2),mem->box_j1->operator[](2));
                    toBB[3]=union_interval(it->box_j1->operator[](3),mem->box_j1->operator[](3));
                    toBB[4]=union_interval(it->box_j1->operator[](4),mem->box_j1->operator[](4));
                   
                    cout << "New (big) toBB is here! " << toBB << endl;
                   
                    return(toBB);
                }
            }
        }
    }
}


//Renvoie les intervalles y0 qui vérifient les conditions de rentrée 
//dans le panier. 
//Paramètre d'entrée: int n: nombre de rebond
//Paramètre d'entrée: IntervalVector init: vecteur de condition initiales
//Paramètre d'entrée: Function f: dynamique
list<Interval> simulate_launch(int n, IntervalVector init, Function f, IntervalVector safe){
    IntervalVector w_in = init; //Conditions initiales amenées à être modifiées 
    IntervalVector toBB = init; //Conditions initiales amenées à être bissectée
    double alpha = 0;           
    double t_cross_exact = 0;   //Temps exact d'intersection entre la balle et le sol
    double duration = 0.0;      //Durée de la simulation
    double Eps = 1e-2;
    double delta_T = 0.0; 
    double c = 0.8;   //elasticité du ballon
    list<Interval> valid_y; //liste des intervalles initiaux valides
    list<IntervalVector> toPlot; //liste des boîtes à afficher
    list<solution_g>::iterator it,mem; //itérateurs sur la liste des solutions et mémoire de la solution précédente
    ivp_ode problem = ivp_ode (f, w_in[4].ub(), w_in); //Problème initial
    simulation simu = simulation (&problem, 0.0, __METH__, __PREC__); //Simulation initiale
    for(int i=0;i<=n;i++){
        //Définition des contracteurs
        //condition y>=0 && vx >= min(vx(t)) &&  x>= min(x(t)) && vy >= min(vy(t))
        //- w_in[2].lb()*all[4] - w_in[0].lb()
        Variable all(5); //x,y,vx,vy,t
        Function A_f_const(all, Return( all[0] - w_in[2].lb()*all[4] - w_in[0].lb()
                                      , all[1] - alpha
                                      , all[2] - w_in[2].lb()
                                      , all[3] - (-9.81*all[4] + w_in[3].lb()) 
                                      , -9.81*pow(all[4],2)/2+w_in[3].ub()*all[4] + w_in[1].ub() -alpha ));  
        NumConstraint A_allConst(A_f_const,GEQ);
        CtcFwdBwd A_ctc_All(A_allConst);
        CtcFixPoint A_fpAll(A_ctc_All, 1e-4);

         //condition y<=eps && vx <= max(vx(t)) &&  x<= max(x(t)) && vy <= max(vy(t))
         // - w_in[2].ub()*all[4] - w_in[0].ub() 
        Function B_f_const(all, Return( all[0] - w_in[2].ub()*all[4] - w_in[0].ub()
                                      , all[1]- alpha
                                      , all[2] - w_in[2].ub() 
                                      , all[3] - (-9.81*all[4] + w_in[3].ub())  
                                      , -9.81*pow(all[4],2)/2+w_in[3].lb()*all[4] + w_in[1].lb() -alpha )); 
        NumConstraint B_allConst(B_f_const,LEQ);
        CtcFwdBwd B_ctc_All(B_allConst);
        CtcFixPoint B_fpAll(B_ctc_All, 1e-7);
        //Définition du problème
        ivp_ode problem = ivp_ode (f, w_in[4].ub(), w_in);
        //Calcul du temps d'intersection exact avec le sol
        t_cross_exact = t_cross_exact + exact_intersect_time(w_in[1].mid(),w_in[3].mid());
        cout << "t_cross_exact: " << t_cross_exact << endl; 
        delta_T = exact_intersect_time(w_in[1].diam(),w_in[3].diam());
        cout << "delta_T: " << delta_T << endl; 
        //MODIFICATION: REMPLACÉ delta_T par 1.
        duration = t_cross_exact+1.;
        cout << "duration of simulation: " << duration << endl; 
        simulation simu = simulation (&problem, duration, __METH__, __PREC__);
        simu.run_simulation(); 
        toBB = get_toBB(&simu, &toPlot, toBB);

        cout << "before first contraction x: " << toBB[0] << endl;
        cout << "before first contraction y : " << toBB[1] << endl;
        cout << "before first contraction vx : " << toBB[2] << endl;
        cout << "before first contraction vy : " << toBB[3] << endl;
        cout << "before first contraction t : " << toBB[4] << endl;
        alpha = 0.;
        A_fpAll.contract(toBB);
        cout << "after first contraction x: " << toBB[0] << endl;
        cout << "after first contraction y : " << toBB[1] << endl;
        cout << "after first contraction vx : " << toBB[2] << endl;
        cout << "after first contraction vy : " << toBB[3] << endl;
        cout << "after first contraction t : " << toBB[4] << endl;
        toPlot.push_back(toBB);
        alpha = Eps;
        B_fpAll.contract(toBB);
        cout << "final x: " << toBB[0] << endl;
        cout << "final y : " << toBB[1] << endl;
        cout << "final vx : " << toBB[2] << endl;
        cout << "final vy : " << toBB[3] << endl;
        cout << "final t : " << toBB[4] << endl;
        toPlot.push_back(toBB);
        w_in[0] = toBB[0];  // x
        w_in[1] = toBB[1];  // y
        w_in[2] = c*toBB[2];  // vx
        w_in[3] = (-c)*toBB[3];  // vy, condition initiale négative pour simuler le rebond
        w_in[4] = toBB[4]; // t
        cout << "Final initial condition: " << w_in << endl;
        //Vérifier les contraintes sur x ici.
        //Si au rebond n-1 on atteint la zone après le panier sans être à la
        //bonne hauteur, on est allé trop loin
        if(is_behind_basket(toBB,safe)){
            std::cerr << "\t\t\033[31mErreur: dépassement du panier\033[0m" << std::endl;
            if (w_in[0].diam() > 0.0001){ //largeur arbitraire
                IntervalVector w_in_1=w_in;
                IntervalVector w_in_2=w_in;
                tuple<Interval,Interval> t = divide_interval(w_in[0]);
                w_in_1[0]=get<0>(t);
                w_in_2[0]=get<1>(t);
                simulate_launch(n,w_in_1,f,safe);
                simulate_launch(n,w_in_2,f,safe);
            }
        }
    }
    //Vérifier les contraintes sur panier ici.
    //if contraintes_y not verified then 
    if( has_crossed_basket(&simu,safe) ){
        cout << "\t\t\033Panier traversé!\033[0m" << endl;
        valid_y.push_back(w_in[1]);
    }
    else{
        if (w_in[1].diam() > 0.0001){ 
            IntervalVector w_in_1=w_in;
            IntervalVector w_in_2=w_in;
            tuple<Interval,Interval> t = divide_interval(w_in[1]);
            w_in_1[1]=get<0>(t);
            w_in_2[1]=get<1>(t);
            simulate_launch(n,w_in_1,f,safe);
            simulate_launch(n,w_in_2,f,safe);
        }
    }
    plot_all(toPlot);
    return(valid_y);
}

int main() {
  
    list<Interval> res;
    //**********
    //conditions initiales
    //**********
    //Vitesse
    Interval Vx0(10.0,10.0);
    Interval Vy0(3.0,3.0);
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
    //Zone du panier 
    IntervalVector basket(5);
    basket[0] = Interval(-0.1, 6.962);
    basket[1] = Interval(3.05,3.05);
    basket[2] = Interval(-400,400); 
    basket[3] = Interval(-400,-0.1); //On ne peut franchir le panier qu'en chutant
    basket[4] = Interval(-0.1,1000); //On peut l'atteindre à n'importe quel temps
    
    //**********
    //dynamique
    //**********
    Variable y(5); 
    Function deriv = Function (y, Return(y[2],y[3],Interval(0.0,0.0),-g,Interval(1.0,1.0)) );
    //*********
    //simulation et affichage des résultats
    //*********
    vibes::beginDrawing ();
    vibes::newFigure("Basket");
    vibes::drawBox(0,6.93,3.05,3.05, "[red]");
    res = simulate_launch(1, State_init,deriv,basket);
    return 0;
}
