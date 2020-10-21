#include <iostream>
#include <math.h>
#include "bindings.h"

using namespace std;

// >>> Main HPolytopeCPP class; compute_volume(), rounding() and generate_samples() volesti methods are included <<<

// Here is the initialization of the HPolytopeCPP class
HPolytopeCPP::HPolytopeCPP() {}
HPolytopeCPP::HPolytopeCPP(double *A_np, double *b_np, int n_hyperplanes, int n_variables){

   MT A;
   VT b;
   A.resize(n_hyperplanes, n_variables);
   b.resize(n_hyperplanes);

   svd_parameters = svd_params(n_variables);

   int index = 0;
   for (int i = 0; i < n_hyperplanes; i++){
      b(i) = b_np[i];
      for (int j=0; j < n_variables; j++){
         A(i,j) = A_np[index];
         index++;
      }
   }

   HP = Hpolytope(n_variables, A, b);

}

// Use a destructor for the HPolytopeCPP object
HPolytopeCPP::~HPolytopeCPP(){}

//////////          Start of "compute_volume"          //////////
double HPolytopeCPP::compute_volume(char* vol_method, char* walk_method, 
                                    int walk_len, double epsilon, int seed){

   double volume;

   if (strcmp(vol_method,"sequence_of_balls") == 0){
      if (strcmp(walk_method,"uniform_ball") == 0){
         volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"CDHR") == 0){
         volume = volume_sequence_of_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"RDHR") == 0){
         volume = volume_sequence_of_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   }
   else if (strcmp(vol_method,"cooling_gaussian") == 0){
      if (strcmp(walk_method,"gaussian_ball") == 0){
         volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_CDHR") == 0){
         volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_RDHR") == 0){
         volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   } else if (strcmp(vol_method,"cooling_balls") == 0){
       if (strcmp(walk_method,"uniform_ball") == 0){
         volume = volume_cooling_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"CDHR") == 0){
         volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"RDHR") == 0){
         volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"billiard") == 0){
         volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, epsilon, walk_len);
       }
   }
   return volume;
}
//////////           End of "compute_volume()"            //////////

//////////         Start of "generate_samples()"          //////////
double HPolytopeCPP::generate_samples(int walk_len, int number_of_points, 
                                      int number_of_points_to_burn, bool boundary,
                                      bool cdhr, bool rdhr, bool gaussian, bool set_L,
                                      bool accelerated_billiard, bool billiard,
                                      bool ball_walk, double a, double L, bool max_ball,
                                      double* inner_point, double radius, double* samples){
   
   RNGType rng(HP.dimension());
   HP.normalize();
   
   int d = HP.dimension();
   Point starting_point; 
   
   // Check for max ball given
   if (max_ball == true){

      VT inner_vec(d);
      for (int i = 0; i < d; i++){
         inner_vec(i) = inner_point[i];
      }
   
      Point inner_point2(inner_vec); 
      CheBall = std::pair<Point, NT>(inner_point2, radius);
      HP.set_InnerBall(CheBall);
      starting_point = inner_point2;
      
   } else {

      //Point default_starting_point = HP.ComputeInnerBall().first;
      starting_point = HP.ComputeInnerBall().first;
   }   
      
   std::list<Point> rand_points;

   if (boundary == true) {
      if (cdhr == true) {
         uniform_sampling_boundary<BCDHRWalk>(rand_points, HP, rng, walk_len,
                                              number_of_points, starting_point,
                                              number_of_points_to_burn);
         } else {
            uniform_sampling_boundary<BRDHRWalk>(rand_points, HP, rng, walk_len, 
                                                 number_of_points, starting_point, 
                                                 number_of_points_to_burn);
         }
   } else if (cdhr == true) {
      if (gaussian == true) {
         gaussian_sampling<GaussianCDHRWalk>(rand_points, HP, rng, walk_len, 
                                             number_of_points, a, starting_point,
                                             number_of_points_to_burn);
      } else {
         uniform_sampling<CDHRWalk>(rand_points, HP, rng, walk_len, number_of_points,
                                    starting_point, number_of_points_to_burn);
      }
   } else if (rdhr == true){
      if (gaussian == true) {
         gaussian_sampling<GaussianRDHRWalk>(rand_points, HP, rng, walk_len, 
                                             number_of_points, a, starting_point, 
                                             number_of_points_to_burn);
      } else {
         uniform_sampling<RDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, 
                                    starting_point, number_of_points_to_burn);
      }
   } else if (billiard == true) {
      if (set_L == true) {
         BilliardWalk WalkType(L);
         uniform_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points,
                          starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<BilliardWalk>(rand_points, HP, rng, walk_len, 
                                        number_of_points, starting_point, 
                                        number_of_points_to_burn);
      }
   } else if (accelerated_billiard == true) {
      if (set_L == true) {
         AcceleratedBilliardWalk WalkType(L);
         uniform_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points,
                          starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<AcceleratedBilliardWalk>(rand_points, HP, rng, walk_len, 
                                        number_of_points, starting_point, 
                                        number_of_points_to_burn);
      }
   } else {
      if (set_L == true) {
         if (gaussian == true) {
            GaussianBallWalk WalkType(L);
            gaussian_sampling(rand_points, HP, rng, WalkType, walk_len,
                              number_of_points, a, starting_point, 
                              number_of_points_to_burn);
            } else {
               BallWalk WalkType(L);
               uniform_sampling(rand_points, HP, rng, WalkType, walk_len,
                                number_of_points, starting_point, 
                                number_of_points_to_burn);
            }
        } else {
            if (gaussian == true) {
               gaussian_sampling<GaussianBallWalk>(rand_points, HP, rng, walk_len, 
                                                   number_of_points, a, starting_point, 
                                                   number_of_points_to_burn);
            } else {
               uniform_sampling<BallWalk>(rand_points, HP, rng, walk_len,
                                          number_of_points, starting_point, 
                                          number_of_points_to_burn);
            }
        }
   }

// The following block of code allows us to parse the matrix with the points we are making
   auto n_si=0;
   for (auto it_s = rand_points.begin(); it_s != rand_points.end(); it_s++){
      for (auto i = 0; i != it_s->dimension(); i++){
         samples[n_si++] = (*it_s)[i];
      }
   }
}
//////////         End of "generate_samples()"          //////////


//////////         Start of "rounding()"          //////////
void HPolytopeCPP::rounding(char* rounding_method, double* new_A, double* new_b,
                            double* T_matrix, double* shift, double &round_value,
                            bool max_ball, double* inner_point, double radius){

   // make a copy of the initial HP which will be used for the rounding step
   auto P(HP);
   RNGType rng(P.dimension());
   P.normalize();
   
   // check for max ball given
   if (max_ball == true ){
      
      // if yes, then read the inner point provided by the user and the radius
      int d = P.dimension();
      VT inner_vec(d);
      
      for (int i = 0; i < d; i++){
         inner_vec(i) = inner_point[i];
      }

      Point inner_point2(inner_vec);
      CheBall = std::pair<Point, NT>(inner_point2, radius);
      P.set_InnerBall(CheBall);

      
   } else if (max_ball == false ) {
      CheBall = P.ComputeInnerBall();
   }
   

   // set the output variable of the rounding step
   round_result round_res;

   // walk length will always be equal to 2
   int walk_len = 2;

   // run the rounding method
   if (strcmp(rounding_method,"min_ellipsoid") == 0){
      round_res = min_sampling_covering_ellipsoid_rounding<AcceleratedBilliardWalk, MT, VT>(P,
                                                                                            CheBall,
                                                                                            walk_len,
                                                                                            rng);
   } else if (strcmp(rounding_method,"svd") == 0){
      round_res = svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, CheBall, walk_len, rng);
   } else if (strcmp(rounding_method, "max_ellipsoid") == 0){
      round_res = max_inscribed_ellipsoid_rounding<MT, VT, NT>(P, CheBall.first);
   }

   // create the new_A matrix
   MT A_to_copy = P.get_mat();
   int n_hyperplanes = P.num_of_hyperplanes();
   int n_variables = P.dimension();

   auto n_si = 0;
   for (int i = 0; i < n_hyperplanes; i++){
      for (int j = 0; j < n_variables; j++){
         new_A[n_si++] = A_to_copy(i,j);
      }
   }

   // create the new_b vector
   VT new_b_temp = P.get_vec();
   for (int i=0; i < n_hyperplanes; i++){
      new_b[i] = new_b_temp[i];
   }

   // create the T matrix
   MT T_matrix_temp = get<0>(round_res);
   auto t_si = 0;
   for (int i = 0; i < n_variables; i++){
      for (int j = 0; j < n_variables; j++){
         T_matrix[t_si++] = T_matrix_temp(i,j);
      }
   }

   // create the shift vector
   VT shift_temp = get<1>(round_res);
   for (int i = 0; i < n_variables; i++){
      shift[i] = shift_temp[i];
   }

   // get the round val value from the rounding step and print int
   round_value = get<2>(round_res);

}
//////////         End of "rounding()"          //////////

double HPolytopeCPP::rounding_svd_step(double* new_A, double* new_b,
                                     double* T_matrix, double* shift,
                                     double* inner_point, double radius) {

   RNGType rng(HP.dimension());
   HP.normalize();
      
   // if yes, then read the inner point provided by the user and the radius
   int d = HP.dimension();
   VT inner_vec(d);
   
   for (int i = 0; i < d; i++){
      inner_vec(i) = inner_point[i];
   }
   
   Point inner_point2(inner_vec);   
   CheBall = std::pair<Point, NT>(inner_point2, radius);   
      
   HP.set_InnerBall(CheBall);   
   unsigned int walk_length = 2;
   svd_rounding_single_step<AcceleratedBilliardWalk, MT, VT>(HP, CheBall, walk_length,
                                                             svd_parameters, rng);

   if (svd_parameters.fail) {
      return -1;
   }

   MT A_to_copy = HP.get_mat();
   int n_hyperplanes = HP.num_of_hyperplanes();
   int n_variables = HP.dimension();

   auto n_si = 0;
   
   for (int i = 0; i < n_hyperplanes; i++){
      for (int j = 0; j < n_variables; j++){
         new_A[n_si++] = A_to_copy(i,j);
      }
   }
   
   // create the new_b vector
   VT new_b_temp = HP.get_vec();
   for (int i=0; i < n_hyperplanes; i++){
      new_b[i] = new_b_temp[i];
   }

   if (svd_parameters.converged) {
      
      // create the T matrix
      MT T_matrix_temp = svd_parameters.T;

      auto t_si = 0;
      for (int i = 0; i < n_variables; i++){
         for (int j = 0; j < n_variables; j++){
            T_matrix[t_si++] = T_matrix_temp(i,j);
         }
      }
   
      // create the shift vector
      VT shift_temp = svd_parameters.T_shift;
      for (int i = 0; i < n_variables; i++){
         shift[i] = shift_temp[i];
      }
      return 1.5;
   }
   return 2.5;
   
}

// >>> The lowDimHPolytopeCPP class; the pre_processing() and the get_full_dimensional_polytope() volesti methods are included <<<

lowDimHPolytopeCPP::lowDimHPolytopeCPP() {}

// Initialize the low dimensional polytope object
lowDimHPolytopeCPP::lowDimHPolytopeCPP(double *A_np, double *b_np, double *A_aeq_np,
                                       double *b_aeq_np, int n_rows_of_A, int n_cols_of_A,
                                       int n_row_of_Aeq, int n_cols_of_Aeq){

   A.resize(n_rows_of_A,n_cols_of_A);
   b.resize(n_rows_of_A);
   
   Aeq.resize(n_row_of_Aeq, n_cols_of_Aeq);
   beq.resize(n_row_of_Aeq);

   int index_1 = 0;
   for (int i = 0; i < n_rows_of_A; i++){
      b(i) = b_np[i];
      for (int j=0; j < n_cols_of_A; j++){
         A(i,j) = A_np[index_1];
         index_1++;
      }
   }

   int index_2 = 0;
   for (int i = 0; i < n_row_of_Aeq; i++){
      beq(i) = b_aeq_np[i];
      for (int j=0; j < n_cols_of_Aeq; j++){
         Aeq(i,j) = A_aeq_np[index_2];
         index_2++;
      }
   }   
}
// Destructor! - never forget about this!
lowDimHPolytopeCPP::~lowDimHPolytopeCPP(){}


// Function to get the full dimensional polytope
int lowDimHPolytopeCPP::full_dimensiolal_polytope(double* N_extra_trans, double* shift,
                                                  double* A_full_extra_trans, double* b_full){
   
   get_full_dim_pol_result result;
   
   // we now run thi C++ function for getting the full dim pol
   result = get_full_dimensional_polytope<Hpolytope>(A, b, Aeq, beq);

   // the outcome of the full_dimensional_polytope()
   Hpolytope full_HP = result.first;
   MT full_HP_A_trans = full_HP.get_mat().transpose();
   VT full_HP_b = full_HP.get_vec();
   MT N_temp_trans = result.second.first.transpose();
   VT shift_temp = result.second.second;   
   
   // Here is what we need to return the output in the Python interface
   // return the full_HP_A matrix to cython
   auto a_si = 0;
   for (int i = 0; i < full_HP_A_trans.rows(); i++){
      for (int j = 0; j < full_HP_A_trans.cols(); j++){
         A_full_extra_trans[a_si++] = full_HP_A_trans(i,j);
      }
   }
   
   // return the full_b matrix to cython
   for (int i=0; i < full_HP_b.rows(); i++){
      b_full[i] = full_HP_b[i];
   }

   // return the N matrix to cython
   auto t_si = 0;
   for (int i = 0; i < N_temp_trans.rows(); i++){
      for (int j = 0; j < N_temp_trans.cols(); j++){
         N_extra_trans[t_si++] = N_temp_trans(i,j);
      }
   }

   // return the shift vector to cython
   for (int i=0; i < shift_temp.rows(); i++){
      shift[i] = shift_temp[i];
   }   
   
   // as we know that N_temp.cols == full_HP_A.cols and likewise for their lines, 
   // we may return just one of those vars
   return N_temp_trans.rows();
} 

