#include <iostream>
#include "HarmonicProblem.h"

using namespace std;

int main()
{
   HarmonicProblem hp;
   hp.ReadFormGrid("data/grid.txt");
   hp.ReadMatrices();
   hp.test = Test(0);

   //for(int i = 0; i < 8; i++)
   //{
   //  cout << i % 2 << " " << (int)floor(i / 2) % 2 << " " << (int)floor(i / 4) % 4 << endl;
   //}

   //for(int i = 0; i < hp.elems_count; i++)
   //{
   //   vector<double> x_nodes_elem(2);       // Координаты конечного элемента по x
   //   vector<double> y_nodes_elem(2);       // Координаты конечного элемента по y
   //   vector<double> z_nodes_elem(2);       // Координаты конечного элемента по y
   //   hp.CalcElemNodes(i, x_nodes_elem, y_nodes_elem, z_nodes_elem);
   //}

   hp.InitializeMemory();
   hp.FormGlobalPortrait();
   hp.AssembleGlobalMatrix();



}