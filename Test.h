#pragma once
#pragma once
using namespace std;

class Test
{
public:

   int N, I, J, L, P;

   Test(const int& n, const int& i, const int& j, const int& l, const int& p) : N(n), I(i), J(j), L(l), P(p) {};

   Test() : I(0), J(0), L(0), P(0) {};

   double us(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return y + z;
         case(2): return y * z;
      };
   }

   double uc(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return 5.0;
         case(1): return x + y;
         case(2): return x * y;
      };
   }

   double fs(const double& x, const double& y, const double& z)
   {
     return -1 * divlambdagradus(x, y, z) +
            -1 * omega() * sigma() * uc(x, y, z) -
            omega() * omega() * chi() * us(x, y, z);
   }

   double fc(const double& x, const double& y, const double& z)
   {
      return -1 * divlambdagraduc(x, y, z) + 
            omega() * sigma() * us(x, y, z) -
            omega() * omega() * chi() * uc(x, y, z);
   }

   double divlambdagradus(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         default: return 0;
      };
   }

   double divlambdagraduc(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         default: return 0;
      };
   }

   double lambda() 
   { 
      switch (I)
      {
         case(0): return 1e2;
         case(1): return 2e3;
         case(2): return 4e4;
         case(3): return 8e5;
      };
   }

   double sigma() 
   { 
      switch (J)
      {
         case(0): return 0;
         case(1): return 1e3;
         case(2): return 1e5;
         case(3): return 1e8;
      };
   }

   double chi() 
   { 
      switch (L)
      {
         case(0): return 8.81e-12;
         case(1): return 1e-12;
         case(2): return 1e-11;
         case(3): return 1e-10;
      };
   }

   double omega() 
   { 
      switch (P)
      {
         case(0): return 1e-4;
         case(1): return 1e1;
         case(2): return 1e5;
         case(3): return 1e9;
      };
   }
};