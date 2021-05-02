#pragma once
#pragma once
using namespace std;

class Test
{
public:

   int N;

   Test(const int& t_N) : N(t_N) {};

   Test() : N(0){};

   double us(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return 2.0;
      };
   }

   double uc(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return 5.0;
      };
   }

   double fs(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return -1 * divgradus(x, y, z) +
            -1 * omega() * sigma() * uc(x, y, z) -
            omega() * omega() * chi() * us(x, y, z);
      };
   }

   double fc(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return -1 * divlambdagraduc(x, y, z) + 
            omega() * sigma() * us(x, y, z) -
            omega() * omega() * chi() * uc(x, y, z);
      };
   }

   double divgradus(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return 0;
      };
   }

   double divlambdagraduc(const double& x, const double& y, const double& z)
   {
      switch(N)
      {
         case(0): return 0;
      };
   }

   double lambda() { return 10; }

   double sigma() { return 1; }

   double chi() { return 1; }

   double omega() { return 2; }
};