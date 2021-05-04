#pragma once
#include <iostream>
#include "Matrix.h"
#include "SLAE.h"
#include "Test.h"
#include <ctime> 

using namespace std;

class HarmonicProblem
{
public:
   Matrix global;                 // ���������� �������
   Matrix fac_global;             // �������� ������������ ���������� �������

   Matrix  GMM, MGM, MMG, MMM;    // ��������������� ������� ��� ���������� ���������
                                  // ��������� ������ � �������� ������ ������

   SLAE slae;                     // �������� ������� ��� ������������������
   SLAE fac_slae;                 // �������� ������� c �������������������

   vector<double> b;              // ���������� ������ ������ �����
   vector<double> loc_b;          // ��������� ������ ������ �����
   vector<double> true_solution;  // ������ �������
   vector<double> solution;       // �������
   vector<int> location;          // ��������� �� ����� ��� ������� ����

   Test test;                     // ���������� � ��������� �������,
                                  // ��������� ������ � �����

   int elems_count;               // ���������� �������� ���������
   int nodes_count;               // ���������� �����

   vector<double> x_nodes;        // ���������� ����� �� x
   vector<double> y_nodes;        // ���������� ����� �� y
   vector<double> z_nodes;        // ���������� ����� �� z
                                  
   int x_nodes_count;             // ���������� ����� �� x
   int y_nodes_count;             // ���������� ����� �� y
   int z_nodes_count;             // ���������� ����� �� z

   vector<vector<int>> regions;   // ���������� � ����������� (��������)
   int regions_count;             // ���������� ��������

   vector<int> x_cords_i;         // ������� �������� ������������ ����� � ������� ����� �� x
   vector<int> y_cords_i;         // ������� �������� ������������ ����� � ������� ����� �� y
   vector<int> z_cords_i;         // ������� �������� ������������ ����� � ������� ����� �� z

   vector<vector<int>> boundaries; // ���������� � ������� ��������
   int bound_count;                // ���������� ������� �������

   // ��������������� ������� ��� ������������ �����
   void ReadFormGridHelp(int& t_nodes_count, vector<double>& t_nodes, vector<int>& t_cords_i, ifstream& fin)
   {
      int t_coords_count;
      fin >> t_coords_count;

      vector<double> t_coords(t_coords_count);

      for(int i = 0; i < t_coords_count; i++)
         fin >> t_coords[i];

      t_nodes_count = 1;
      t_nodes = vector<double>(1);
      t_nodes[0] = t_coords[0];
      vector<int> n_t(t_coords_count - 1);

      for(int i = 0; i < t_coords_count - 1; i++)
      {
         fin >> n_t[i];
         t_nodes.resize(t_nodes_count + n_t[i]);

         double h = (t_coords[i + 1] - t_coords[i]) / n_t[i];

         for(int j = 0; j < n_t[i]; j++)
            t_nodes[j + t_nodes_count] = t_nodes[t_nodes_count - 1] + h * (j + 1);

         t_nodes_count += n_t[i];
      }

      // �������� �������� ������ �����������
      for(int i = 1; i < n_t.size(); i++)
         n_t[i] += n_t[i - 1];

      t_cords_i = vector<int>(n_t.size() + 1);

      for(int i = 0; i < n_t.size(); i++)
         t_cords_i[i + 1] = n_t[i];

      elems_count *= n_t[n_t.size() - 1];
   }

   // ���������� � ������������ ����� �� ����� file_name
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);
      string fake;

      elems_count = 1;

      // ���������� ������������ ����� � ������������ ����� �� x
      fin >> fake;
      ReadFormGridHelp(x_nodes_count, x_nodes, x_cords_i, fin);

      // ���������� ������������ ����� � ������������ ����� �� y
      fin >> fake;
      ReadFormGridHelp(y_nodes_count, y_nodes, y_cords_i, fin);

      // ���������� ������������ ����� � ������������ ����� �� z
      fin >> fake;
      ReadFormGridHelp(z_nodes_count, z_nodes, z_cords_i, fin);
      
      // ���������� ���������� � �����������
      fin >> fake;
      fin >> regions_count;
      regions = vector<vector<int>>(regions_count, vector<int>(6));

      for(int i = 0; i < regions_count; i++)
         fin >> regions[i][0] >> regions[i][1] >> regions[i][2] >> regions[i][3] >> regions[i][4] >> regions[i][5];

      fin.close();

      // ���������� �������� ������ ����������� � ������������ � ���������� �����
      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         regions[reg_i][0] = x_cords_i[regions[reg_i][0]];
         regions[reg_i][1] = x_cords_i[regions[reg_i][1]];
         regions[reg_i][2] = y_cords_i[regions[reg_i][2]];
         regions[reg_i][3] = y_cords_i[regions[reg_i][3]];
         regions[reg_i][4] = z_cords_i[regions[reg_i][4]];
         regions[reg_i][5] = z_cords_i[regions[reg_i][5]];
      }

      nodes_count = x_nodes_count * y_nodes_count * z_nodes_count;
   }

   // ���������� ���������� � ������� �������� �� ����� file_name
   void ReadBoundaries(const string& file_name)
   {
      ifstream fin(file_name);

      fin >> bound_count;

      boundaries = vector<vector<int>>(bound_count, vector<int>(7));

      for(int bound_i = 0; bound_i < bound_count; bound_i++)
      {
         for(int i = 0; i < 7; i++)
            fin >> boundaries[bound_i][i];

         boundaries[bound_i][1] = x_cords_i[boundaries[bound_i][1]];
         boundaries[bound_i][2] = x_cords_i[boundaries[bound_i][2]];
         boundaries[bound_i][3] = y_cords_i[boundaries[bound_i][3]];
         boundaries[bound_i][4] = y_cords_i[boundaries[bound_i][4]];
         boundaries[bound_i][5] = z_cords_i[boundaries[bound_i][5]];
         boundaries[bound_i][6] = z_cords_i[boundaries[bound_i][6]];
      }

      fin.close();
   }

   // ���������� ��������������� ������ ��� ������������
   // ������ ��������� � �����
   void ReadMatrices()
   {
      GMM = Matrix(8);
      GMM.ReadDiTr("data/GMM.txt");

      MGM = Matrix(8);
      MGM.ReadDiTr("data/MGM.txt");

      MMG = Matrix(8);
      MMG.ReadDiTr("data/MMG.txt");

      MMM = Matrix(8);
      MMM.ReadDiTr("data/MMM.txt");
   }

   // ��������� ������ ��� �������
   void InitializeMemory()
   {
      slae = SLAE(nodes_count * 2, 10000, 1e-20);
      fac_slae = SLAE(nodes_count * 2, 10000, 1e-20);

      global.ind = vector<int>(nodes_count * 2 + 1);
      b = vector<double>(nodes_count * 2);

      solution = vector<double>(nodes_count * 2);
      true_solution = vector<double>(nodes_count * 2);
      location = vector<int>(nodes_count * 2);

      fac_global = Matrix(nodes_count * 2, 0);
   }

   // ���������� ������� global_indices ���������, ���������������� ���������� ���������
   // ����� ��������� �������� � ������� elem_index(���������� � ����)
   void CalcGlobalIndices(int elem_i, vector<int>& global_indices)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int k = elem_i % (x_nodes_count - 1) + x_nodes_count * floor(elem_i / (x_nodes_count - 1));
      k = k % (x_nodes_count * (y_nodes_count - 1)) + (x_nodes_count * y_nodes_count) * floor(k / (x_nodes_count * (y_nodes_count - 1)));

      global_indices[0] = (k + 0);
      global_indices[1] = (k + 1);

      global_indices[2] = (k + x_nodes_count + 0);
      global_indices[3] = (k + x_nodes_count + 1);

      global_indices[4] = (k + x_nodes_count * y_nodes_count + 0);
      global_indices[5] = (k + x_nodes_count * y_nodes_count + 1);

      global_indices[6] = (k + x_nodes_count * y_nodes_count + x_nodes_count + 0);
      global_indices[7] = (k + x_nodes_count * y_nodes_count + x_nodes_count + 1);
   }

   // ����� ������� �� ������ ��������� ��������
   int CalcRegionIndex(const int& elem_i)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int x0 = (elem_i) % (n_coords - 1) * 2 + 1;
      int y0 = floor((elem_i) / (n_coords - 1)) * 2 + 1;

      int found_reg_i = -1;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         if(regions[reg_i][0] <= x0 && x0 <= regions[reg_i][1] &&
            regions[reg_i][2] <= y0 && y0 <= regions[reg_i][3])
         {
            found_reg_i = reg_i;
            break;
         }
      }
      return found_reg_i;
   }

   // ��������������� ������� ��� ������������ ��������
   void IncertToRow(Matrix& mat, const int& r, const int& c)
   {
      int i_in_jg = mat.ind[r];
      int prof_len = mat.ind[r + 1] - mat.ind[r];

      bool found = false;

      for(int k = i_in_jg; k < i_in_jg + prof_len; k++)
         if(mat.columns_ind[k] == c)
         {
            found = true;
            break;
         }

      if(!found)
      {
         for(int l = r + 1; l < mat.ind.size(); l++)
            mat.ind[l]++;

         int k = i_in_jg;

         while((k < i_in_jg + prof_len) && mat.columns_ind[k] < c)
            k++;

         mat.columns_ind.insert(mat.columns_ind.begin() + k, c);
      }
   }

   // ������������ �������� ���������� ������� � ������ ������� ���������
   void FormGlobalPortrait()
   {
      global.ind[0] = global.ind[1] = 0;

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         int reg_i = CalcRegionIndex(elem_i);

         //if(reg_i != -1)
         {
            vector<int> global_indices(8);
            CalcGlobalIndices(elem_i, global_indices);
            vector<vector<vector<int>>> help(8);

            for(int i = 0; i < 8; i++)
               help[i].resize(8);

            for(int i = 0; i < 8; i++)
               for(int j = 0; j < 8; j++)
                  help[i][j] = { global_indices[i], global_indices[j] };

            for(int i = 0; i < 8; i++)
            {
               IncertToRow(global, help[i][i][0] * 2 + 1, help[i][i][1] * 2);
            }

            for(int i = 0; i < 8; i++)
               for(int j = 0; j < i; j++)
               {
                  IncertToRow(global, help[i][j][0] * 2, help[i][j][1] * 2);
                  IncertToRow(global, help[i][j][0] * 2, help[i][j][1] * 2 + 1);
                  IncertToRow(global, help[i][j][0] * 2 + 1, help[i][j][1] * 2);
                  IncertToRow(global, help[i][j][0] * 2 + 1, help[i][j][1] * 2 + 1);
               }
         }
      }

      global.size = global.ind.size() - 1;
      global.diag.resize(global.size);
      
      global.tr_size = global.columns_ind.size();
      global.bot_tr.resize(global.tr_size);
      global.top_tr.resize(global.tr_size);
   }

   // ������������ �������� ������� mat
   void FormPortrait(Matrix& mat)
   {
      mat.ind[0] = mat.ind[1] = 0;

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         int reg_i = CalcRegionIndex(elem_i);

         //if(reg_i != -1)
         {
            vector<int> mat_indices(8);
            CalcGlobalIndices(elem_i, mat_indices);
            vector<vector<vector<int>>> help(8);

            for(int i = 0; i < 8; i++)
               help[i].resize(8);

            for(int i = 0; i < 8; i++)
               for(int j = 0; j < 8; j++)
                  help[i][j] = { mat_indices[i], mat_indices[j] };

            for(int i = 0; i < 8; i++)
               for(int j = 0; j < i; j++)
                  IncertToRow(mat, help[i][j][0], help[i][j][1]);
         }
      }

      mat.size = mat.ind.size() - 1;
      mat.diag.resize(mat.size);

      mat.tr_size = mat.columns_ind.size();
      mat.bot_tr.resize(mat.tr_size);
      mat.top_tr.resize(mat.tr_size);
   }

   // ���������� � ������� ������� � �������� [row_i][col_i] �������� val
   void AddToMat(Matrix& mat, const int& row_i, const int& col_i, const double& val_l, const double& val_u)
   {
      if(row_i == col_i)
         mat.diag[row_i] += val_l;
      else
      {
         int beg_prof = mat.ind[row_i];
         int end_prof = mat.ind[row_i + 1];

         for(int i_in_prof = beg_prof; i_in_prof < end_prof; i_in_prof++)
         {
            if(mat.columns_ind[i_in_prof] == col_i)
            {
               mat.bot_tr[i_in_prof] += val_l;
               mat.top_tr[i_in_prof] += val_u;
               break;
            }
         }
      }
   }

   // ������� ���������� ����� ��������� ��������
   // � ������� elem_index(���������� � ����)
   void CalcElemNodes(const int& elem_i, vector<double>& x_elem_nodes, vector<double>& y_elem_nodes, vector<double>& z_elem_nodes)
   {
      int n_x = x_nodes_count - 1;
      int n_y = y_nodes_count - 1;
      int n_z = z_nodes_count - 1;

      int x0 = elem_i % n_x;
      int y0 = (int)floor(elem_i / n_x) % n_y;
      int z0 = (int)floor(elem_i / (n_x * n_y));

      x_elem_nodes[0] = x_nodes[x0];
      x_elem_nodes[1] = x_nodes[x0 + 1];

      y_elem_nodes[0] = y_nodes[y0];
      y_elem_nodes[1] = y_nodes[y0 + 1];

      z_elem_nodes[0] = z_nodes[z0];
      z_elem_nodes[1] = z_nodes[z0 + 1];

      //cout << x0 << " " << y0 << " " << z0 << endl;
   }

   // ������ ���������� �������
   void AssembleGlobalMatrix()
   {
      vector<int> global_indices(8);

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         int reg_i = CalcRegionIndex(elem_i);
         CalcGlobalIndices(elem_i, global_indices);

         //if(reg_i == -1)
         //{
         //   for(int i = 0; i < 8; i++)
         //   {
         //      global.diag[global_indices[i] * 2    ] = 1;
         //      global.diag[global_indices[i] * 2 + 1] = 1;

         //      //sigma_mass_mat.diag[global_indices[i]] = 0;
         //      //sigma_mass_mat.diag[global_indices[i] + 1] = 0;

         //      //chi_mass_mat.diag[global_indices[i]] = 0;
         //      //chi_mass_mat.diag[global_indices[i] + 1] = 0;

         //      b[global_indices[i] * 2    ] = 0;
         //      b[global_indices[i] * 2 + 1] = 0;

         //      location[global_indices[i] * 2    ] = 2;
         //      location[global_indices[i] * 2 + 1] = 2;
         //   }
         //}
         //else
         {
            vector<double> x_elem_nodes(2);       // ���������� ��������� �������� �� x
            vector<double> y_elem_nodes(2);       // ���������� ��������� �������� �� y
            vector<double> z_elem_nodes(2);       // ���������� ��������� �������� �� z
            CalcElemNodes(elem_i, x_elem_nodes, y_elem_nodes, z_elem_nodes);

            double hx = x_elem_nodes[1] - x_elem_nodes[0];
            double hy = y_elem_nodes[1] - y_elem_nodes[0];
            double hz = z_elem_nodes[1] - z_elem_nodes[0];

            vector<double> local_fs(8);
            vector<double> local_fc(8);

            for(int i = 0; i < 8; i++)
            {
               double x = x_elem_nodes[i % 2];
               double y = y_elem_nodes[(int)floor(i / 2) % 2];
               double z = z_elem_nodes[(int)floor(i / 4) % 4];

               local_fs[i] = test.fs(x, y, z);
               local_fc[i] = test.fc(x, y, z);

               true_solution[global_indices[i] * 2    ] = test.us(x, y, z);
               true_solution[global_indices[i] * 2 + 1] = test.uc(x, y, z);

               for(int j = 0; j <= i; j++)
               {
                  double p = 0;
                  double c = 0;

                  if(i == j)
                  {
                     p = test.lambda() * (hy * hz / (hx * 36) * GMM.diag[i] + hx * hz / (hy * 36) * MGM.diag[i] + hx * hy / (hz * 36) * MMG.diag[i]);
                     p -= test.omega() * test.omega() * test.chi() * hx * hy * hz / 216.0 * MMM.diag[i];
                     c = test.omega() * test.sigma() * hx * hy * hz / 216.0 * MMM.diag[i];
                  }
                  else
                  {
                     int tr_i = MMM.ind[i] + j;

                     p = test.lambda() * (hy * hz / (hx * 36) * GMM.bot_tr[tr_i] + hx * hz / (hy * 36) * MGM.bot_tr[tr_i] + hx * hy / (hz * 36) * MMG.bot_tr[tr_i]);
                     p -= test.omega() * test.omega() * test.chi() * hx * hy * hz / 216.0 * MMM.bot_tr[tr_i];
                     c = test.omega() * test.sigma() * hx * hy * hz / 216.0 * MMM.bot_tr[tr_i];
                  }

                  AddToMat(global, global_indices[i] * 2    , global_indices[j] * 2    , p, p);
                  AddToMat(global, global_indices[i] * 2 + 1, global_indices[j] * 2 + 1, p, p);

                  AddToMat(global, global_indices[i] * 2 + 1, global_indices[j] * 2    , c, -c);
                  AddToMat(global, global_indices[i] * 2    , global_indices[j] * 2 + 1, -c, c);
               }
            }

            vector<double> local_bs(8);
            vector<double> local_bc(8);
            MMM.MatVecMult(local_fs, local_bs, MMM.bot_tr, MMM.top_tr);
            MMM.MatVecMult(local_fc, local_bc, MMM.bot_tr, MMM.top_tr);

            for(int i = 0; i < 8; i++)
            {
               b[global_indices[i] * 2    ] += hx * hy * hz / 216.0 * local_bs[i];
               b[global_indices[i] * 2 + 1] += hx * hy * hz / 216.0 * local_bc[i];
            }
         }
      }
   }

   // ���� ������ ������� ������� �� ������� � ������� line_i
   void FirstBoundOnLine(int line_i, const double& x, const double& y, const double& z)
   {
      line_i *= 2;

      global.diag[line_i    ] = 1;
      global.diag[line_i + 1] = 1;

      true_solution[line_i    ] = test.us(x, y, z);
      true_solution[line_i + 1] = test.uc(x, y, z);

      b[line_i    ] = test.us(x, y, z);
      b[line_i + 1] = test.uc(x, y, z);

      for(int prof_i = global.ind[line_i]; prof_i < global.ind[line_i + 1]; prof_i++)
         global.bot_tr[prof_i] = 0;

      for(int prof_i = global.ind[line_i + 1]; prof_i < global.ind[line_i + 2]; prof_i++)
         global.bot_tr[prof_i] = 0;

      for(int i = 0; i < nodes_count * 2; i++)
      {
         for(int prof_i = global.ind[i]; prof_i < global.ind[i + 1]; prof_i++)
         {
            if(global.columns_ind[prof_i] == line_i || global.columns_ind[prof_i] == line_i + 1)
               global.top_tr[prof_i] = 0;
         }
      }
   }

   // ���� ������ ������� �������
   void AccountFirstBound()
   {
      for(int x_i = 0; x_i < x_nodes_count; x_i++)
      {
         for(int y_i = 0; y_i < y_nodes_count; y_i++)
         {
            for(int z_i = 0; z_i < z_nodes_count; z_i++)
            {
               for(int bound_i = 0; bound_i < bound_count; bound_i++)
               {
                  if(boundaries[bound_i][1] <= x_i && x_i <= boundaries[bound_i][2] &&
                     boundaries[bound_i][3] <= y_i && y_i <= boundaries[bound_i][4] &&
                     boundaries[bound_i][5] <= z_i && z_i <= boundaries[bound_i][6])
                  {
                     int i = x_i + y_i * x_nodes_count + z_i * x_nodes_count * y_nodes_count;
                     FirstBoundOnLine(i, x_nodes[x_i], y_nodes[y_i], z_nodes[z_i]);
                     location[i * 2   ] = 1;
                     location[i * 2 + 1] = 1;
                     break;
                  }
               }
            }
         }
      }
   }

   // ���������� �������
   void Solve(ofstream& fout)
   {
      slae.b = b;
      global.DiagFact(fac_global);

      vector<double> x0(nodes_count * 2);

      double start_time = clock();
      int iter = slae.ConjGradPredMethod(x0, solution, global, fac_slae, fac_global);
      double end_time = clock();

      fout << iter << "\t" << (end_time - start_time) / CLOCKS_PER_SEC * 1000 << "\t";
   }

   // ����� ������� �� ��������� ���� t � ����� fout 
   void PrintSolution(ofstream& fout)
   {
      /*fout << setw(14) << "x" << setw(14) << "y" << setw(14) << "z";
      fout << setw(14) << "prec" << setw(14) << "calc" << setw(14) << "diff";
      fout << setw(5) << "n" << " loc" << endl;*/

      double norm = 0, norm_u = 0;
      for(int z_i = 0; z_i < z_nodes_count; z_i++)
      {
         for(int y_i = 0; y_i < y_nodes_count; y_i++)
         {
            for(int x_i = 0; x_i < x_nodes_count; x_i++)
            {
               int i = x_i + y_i * x_nodes_count + z_i * x_nodes_count * y_nodes_count;

               double prec = true_solution[i * 2];
               double calc = solution[i * 2];

               //if(x_i % 32 == 0 && y_i % 32 == 0)
              /* {
                  fout << scientific;
                  fout << setw(14) << x_nodes[x_i];
                  fout << setw(14) << y_nodes[y_i];
                  fout << setw(14) << z_nodes[z_i];

                  fout << setw(14) << prec;
                  fout << setw(14) << calc;
                  fout << setw(14) << abs(true_solution[i * 2] - solution[i * 2]);

                  fout << fixed << setw(5) << i * 2;

                  if(location[i * 2] == 2)
                     fout << " outer";
                  else if(location[i * 2] == 1)
                     fout << " border";
                  else
                     fout << " inner";

                  fout << endl;

               }*/

               prec = true_solution[i * 2 + 1];
               calc = solution[i * 2 + 1];

               //if(x_i % 32 == 0 && y_i % 32 == 0)
               /*{
                  fout << scientific;
                  fout << setw(14) << x_nodes[x_i];
                  fout << setw(14) << y_nodes[y_i];
                  fout << setw(14) << z_nodes[z_i];

                  fout << setw(14) << prec;
                  fout << setw(14) << calc;
                  fout << setw(14) << abs(true_solution[i * 2 + 1] - solution[i * 2 + 1]);

                  fout << fixed << setw(5) << i * 2 + 1;

                  if(location[i * 2 + 1] == 2)
                     fout << " outer";
                  else if(location[i * 2 + 1] == 1)
                     fout << " border";
                  else
                     fout << " inner";

                  fout << endl;

               }*/
               norm_u += prec * prec;
               norm += abs(prec - calc) * abs(prec - calc);
            }
         }
      }

      /*fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*|| = " << scientific << sqrt(norm) << endl;*/
      fout << scientific << sqrt(norm) << "\t";
      fout << scientific << sqrt(norm) / sqrt(norm_u) << endl;
   }
};