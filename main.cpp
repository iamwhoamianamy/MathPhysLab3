#include <iostream>
#include "HarmonicProblem.h"

using namespace std;

int main()
{
   HarmonicProblem hp;
   hp.ReadFormGrid("data/grid.txt");

   for(int i = 0; i < hp.elems_count; i++)
   {
      vector<int> vec(8);
      hp.CalcGlobalIndices(i, vec);
   }

}