#include <iostream>
#include "HarmonicProblem.h"

using namespace std;

int main()
{
   HarmonicProblem hp;
   hp.ReadFormGrid("data/grid.txt");
   hp.ReadBoundaries("data/boundaries.txt");
   hp.ReadMatrices();
   hp.InitializeMemory();
   hp.FormGlobalPortrait();

   ofstream fout("result.txt");

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         for (int l = 0; l < 4; l++)
            for (int p = 0; p < 4; p++)
            {
               hp.test = Test(2, i, j, l, p);
               hp.AssembleGlobalMatrix();
               hp.AccountFirstBound();
               hp.Solve(fout);
               fout << p << "\t" << i << "\t" << j << "\t" << l << "\t";
            }

   hp.PrintSolution(fout);
   fout.close();
}