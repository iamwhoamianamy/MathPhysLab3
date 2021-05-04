#include <iostream>
#include "HarmonicProblem.h"

using namespace std;

int main()
{
   HarmonicProblem hp;
   hp.ReadFormGrid("data/grid.txt");
   hp.ReadBoundaries("data/boundaries.txt");
   hp.ReadMatrices();
   //hp.InitializeMemory();
   //hp.FormGlobalPortrait();

   ofstream fout("result.txt");

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         for (int l = 0; l < 4; l++)
            for (int p = 0; p < 4; p++)
            {
               hp.InitializeMemory();
               hp.FormGlobalPortrait();
               hp.test = Test(1, i, j, l, p);
               hp.AssembleGlobalMatrix();
               hp.AccountFirstBound();
               fout << p << "\t" << i << "\t" << j << "\t" << l << "\t" << endl;
               hp.Solve(fout);
               fout << endl;
               hp.PrintSolution(fout);
               cout << p << "\t" << i << "\t" << j << "\t" << l << "\t" << endl;
            }

   //hp.test = Test(1, 1, 0, 0, 0);
   //hp.AssembleGlobalMatrix();
   //hp.AccountFirstBound();
   ////fout << p << "\t" << i << "\t" << j << "\t" << l << "\t" << endl;
   //hp.Solve(fout);
   //fout << endl;
   //hp.PrintSolution(fout);
   ////cout << p << "\t" << i << "\t" << j << "\t" << l << "\t" << endl;


   fout.close();
}