#include <iostream>
#include "HarmonicProblem.h"

using namespace std;

int main()
{
   HarmonicProblem hp;
   hp.ReadFormGrid("data/grid.txt");
   hp.ReadBoundaries("data/boundaries.txt");
   hp.ReadMatrices();

   ofstream fout("result.txt");

   //for (int p = 0; p < 3; p++)
      //for (int i = 0; i < 3; i++)
         //for (int j = 0; j < 3; j++)
            //for (int l = 0; l < 3; l++)
            {
               hp.InitializeMemory();
               hp.FormGlobalPortrait();
              // hp.test = Test(1, i, j, l, p);
               hp.test = Test(1, 0, 2, 0, 2);
               hp.AssembleGlobalMatrix();
               hp.AccountFirstBound();
               hp.test.Print(fout);
               hp.Solve(fout);
               hp.PrintSolution(fout);
            }

   fout.close();
}