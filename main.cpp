#include <iostream>
#include "HarmonicProblem.h"

using namespace std;

int main()
{
   HarmonicProblem hp;
   hp.ReadFormGrid("data/grid.txt");
   hp.ReadBoundaries("data/boundaries.txt");
   hp.ReadMatrices();

   hp.test = Test(0);
   hp.InitializeMemory();

   hp.FormGlobalPortrait();
   hp.AssembleGlobalMatrix();
   hp.AccountFirstBound();
   hp.Solve();

   ofstream fout("result.txt");

   hp.PrintSolution(fout);

   fout.close();


}