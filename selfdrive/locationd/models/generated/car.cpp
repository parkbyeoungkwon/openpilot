
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6439155458356320514) {
   out_6439155458356320514[0] = delta_x[0] + nom_x[0];
   out_6439155458356320514[1] = delta_x[1] + nom_x[1];
   out_6439155458356320514[2] = delta_x[2] + nom_x[2];
   out_6439155458356320514[3] = delta_x[3] + nom_x[3];
   out_6439155458356320514[4] = delta_x[4] + nom_x[4];
   out_6439155458356320514[5] = delta_x[5] + nom_x[5];
   out_6439155458356320514[6] = delta_x[6] + nom_x[6];
   out_6439155458356320514[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3744964575875050616) {
   out_3744964575875050616[0] = -nom_x[0] + true_x[0];
   out_3744964575875050616[1] = -nom_x[1] + true_x[1];
   out_3744964575875050616[2] = -nom_x[2] + true_x[2];
   out_3744964575875050616[3] = -nom_x[3] + true_x[3];
   out_3744964575875050616[4] = -nom_x[4] + true_x[4];
   out_3744964575875050616[5] = -nom_x[5] + true_x[5];
   out_3744964575875050616[6] = -nom_x[6] + true_x[6];
   out_3744964575875050616[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_6996129398019384578) {
   out_6996129398019384578[0] = 1.0;
   out_6996129398019384578[1] = 0.0;
   out_6996129398019384578[2] = 0.0;
   out_6996129398019384578[3] = 0.0;
   out_6996129398019384578[4] = 0.0;
   out_6996129398019384578[5] = 0.0;
   out_6996129398019384578[6] = 0.0;
   out_6996129398019384578[7] = 0.0;
   out_6996129398019384578[8] = 0.0;
   out_6996129398019384578[9] = 1.0;
   out_6996129398019384578[10] = 0.0;
   out_6996129398019384578[11] = 0.0;
   out_6996129398019384578[12] = 0.0;
   out_6996129398019384578[13] = 0.0;
   out_6996129398019384578[14] = 0.0;
   out_6996129398019384578[15] = 0.0;
   out_6996129398019384578[16] = 0.0;
   out_6996129398019384578[17] = 0.0;
   out_6996129398019384578[18] = 1.0;
   out_6996129398019384578[19] = 0.0;
   out_6996129398019384578[20] = 0.0;
   out_6996129398019384578[21] = 0.0;
   out_6996129398019384578[22] = 0.0;
   out_6996129398019384578[23] = 0.0;
   out_6996129398019384578[24] = 0.0;
   out_6996129398019384578[25] = 0.0;
   out_6996129398019384578[26] = 0.0;
   out_6996129398019384578[27] = 1.0;
   out_6996129398019384578[28] = 0.0;
   out_6996129398019384578[29] = 0.0;
   out_6996129398019384578[30] = 0.0;
   out_6996129398019384578[31] = 0.0;
   out_6996129398019384578[32] = 0.0;
   out_6996129398019384578[33] = 0.0;
   out_6996129398019384578[34] = 0.0;
   out_6996129398019384578[35] = 0.0;
   out_6996129398019384578[36] = 1.0;
   out_6996129398019384578[37] = 0.0;
   out_6996129398019384578[38] = 0.0;
   out_6996129398019384578[39] = 0.0;
   out_6996129398019384578[40] = 0.0;
   out_6996129398019384578[41] = 0.0;
   out_6996129398019384578[42] = 0.0;
   out_6996129398019384578[43] = 0.0;
   out_6996129398019384578[44] = 0.0;
   out_6996129398019384578[45] = 1.0;
   out_6996129398019384578[46] = 0.0;
   out_6996129398019384578[47] = 0.0;
   out_6996129398019384578[48] = 0.0;
   out_6996129398019384578[49] = 0.0;
   out_6996129398019384578[50] = 0.0;
   out_6996129398019384578[51] = 0.0;
   out_6996129398019384578[52] = 0.0;
   out_6996129398019384578[53] = 0.0;
   out_6996129398019384578[54] = 1.0;
   out_6996129398019384578[55] = 0.0;
   out_6996129398019384578[56] = 0.0;
   out_6996129398019384578[57] = 0.0;
   out_6996129398019384578[58] = 0.0;
   out_6996129398019384578[59] = 0.0;
   out_6996129398019384578[60] = 0.0;
   out_6996129398019384578[61] = 0.0;
   out_6996129398019384578[62] = 0.0;
   out_6996129398019384578[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8382255420543571790) {
   out_8382255420543571790[0] = state[0];
   out_8382255420543571790[1] = state[1];
   out_8382255420543571790[2] = state[2];
   out_8382255420543571790[3] = state[3];
   out_8382255420543571790[4] = state[4];
   out_8382255420543571790[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8382255420543571790[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8382255420543571790[7] = state[7];
}
void F_fun(double *state, double dt, double *out_3259366650936592814) {
   out_3259366650936592814[0] = 1;
   out_3259366650936592814[1] = 0;
   out_3259366650936592814[2] = 0;
   out_3259366650936592814[3] = 0;
   out_3259366650936592814[4] = 0;
   out_3259366650936592814[5] = 0;
   out_3259366650936592814[6] = 0;
   out_3259366650936592814[7] = 0;
   out_3259366650936592814[8] = 0;
   out_3259366650936592814[9] = 1;
   out_3259366650936592814[10] = 0;
   out_3259366650936592814[11] = 0;
   out_3259366650936592814[12] = 0;
   out_3259366650936592814[13] = 0;
   out_3259366650936592814[14] = 0;
   out_3259366650936592814[15] = 0;
   out_3259366650936592814[16] = 0;
   out_3259366650936592814[17] = 0;
   out_3259366650936592814[18] = 1;
   out_3259366650936592814[19] = 0;
   out_3259366650936592814[20] = 0;
   out_3259366650936592814[21] = 0;
   out_3259366650936592814[22] = 0;
   out_3259366650936592814[23] = 0;
   out_3259366650936592814[24] = 0;
   out_3259366650936592814[25] = 0;
   out_3259366650936592814[26] = 0;
   out_3259366650936592814[27] = 1;
   out_3259366650936592814[28] = 0;
   out_3259366650936592814[29] = 0;
   out_3259366650936592814[30] = 0;
   out_3259366650936592814[31] = 0;
   out_3259366650936592814[32] = 0;
   out_3259366650936592814[33] = 0;
   out_3259366650936592814[34] = 0;
   out_3259366650936592814[35] = 0;
   out_3259366650936592814[36] = 1;
   out_3259366650936592814[37] = 0;
   out_3259366650936592814[38] = 0;
   out_3259366650936592814[39] = 0;
   out_3259366650936592814[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3259366650936592814[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3259366650936592814[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3259366650936592814[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3259366650936592814[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3259366650936592814[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3259366650936592814[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3259366650936592814[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3259366650936592814[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3259366650936592814[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3259366650936592814[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3259366650936592814[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3259366650936592814[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3259366650936592814[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3259366650936592814[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3259366650936592814[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3259366650936592814[56] = 0;
   out_3259366650936592814[57] = 0;
   out_3259366650936592814[58] = 0;
   out_3259366650936592814[59] = 0;
   out_3259366650936592814[60] = 0;
   out_3259366650936592814[61] = 0;
   out_3259366650936592814[62] = 0;
   out_3259366650936592814[63] = 1;
}
void h_25(double *state, double *unused, double *out_2795798451157236819) {
   out_2795798451157236819[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3511382771814294907) {
   out_3511382771814294907[0] = 0;
   out_3511382771814294907[1] = 0;
   out_3511382771814294907[2] = 0;
   out_3511382771814294907[3] = 0;
   out_3511382771814294907[4] = 0;
   out_3511382771814294907[5] = 0;
   out_3511382771814294907[6] = 1;
   out_3511382771814294907[7] = 0;
}
void h_24(double *state, double *unused, double *out_3100482985859878334) {
   out_3100482985859878334[0] = state[4];
   out_3100482985859878334[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7979568905426377477) {
   out_7979568905426377477[0] = 0;
   out_7979568905426377477[1] = 0;
   out_7979568905426377477[2] = 0;
   out_7979568905426377477[3] = 0;
   out_7979568905426377477[4] = 1;
   out_7979568905426377477[5] = 0;
   out_7979568905426377477[6] = 0;
   out_7979568905426377477[7] = 0;
   out_7979568905426377477[8] = 0;
   out_7979568905426377477[9] = 0;
   out_7979568905426377477[10] = 0;
   out_7979568905426377477[11] = 0;
   out_7979568905426377477[12] = 0;
   out_7979568905426377477[13] = 1;
   out_7979568905426377477[14] = 0;
   out_7979568905426377477[15] = 0;
}
void h_30(double *state, double *unused, double *out_7887680365731642356) {
   out_7887680365731642356[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5321462390946073429) {
   out_5321462390946073429[0] = 0;
   out_5321462390946073429[1] = 0;
   out_5321462390946073429[2] = 0;
   out_5321462390946073429[3] = 0;
   out_5321462390946073429[4] = 1;
   out_5321462390946073429[5] = 0;
   out_5321462390946073429[6] = 0;
   out_5321462390946073429[7] = 0;
}
void h_26(double *state, double *unused, double *out_3123920273549208447) {
   out_3123920273549208447[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2193210383124914366) {
   out_2193210383124914366[0] = 0;
   out_2193210383124914366[1] = 0;
   out_2193210383124914366[2] = 0;
   out_2193210383124914366[3] = 0;
   out_2193210383124914366[4] = 0;
   out_2193210383124914366[5] = 0;
   out_2193210383124914366[6] = 0;
   out_2193210383124914366[7] = 1;
}
void h_27(double *state, double *unused, double *out_8394801905458863088) {
   out_8394801905458863088[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4033880403109448117) {
   out_4033880403109448117[0] = 0;
   out_4033880403109448117[1] = 0;
   out_4033880403109448117[2] = 0;
   out_4033880403109448117[3] = 1;
   out_4033880403109448117[4] = 0;
   out_4033880403109448117[5] = 0;
   out_4033880403109448117[6] = 0;
   out_4033880403109448117[7] = 0;
}
void h_29(double *state, double *unused, double *out_5557390440317063615) {
   out_5557390440317063615[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1497663621817241476) {
   out_1497663621817241476[0] = 0;
   out_1497663621817241476[1] = 1;
   out_1497663621817241476[2] = 0;
   out_1497663621817241476[3] = 0;
   out_1497663621817241476[4] = 0;
   out_1497663621817241476[5] = 0;
   out_1497663621817241476[6] = 0;
   out_1497663621817241476[7] = 0;
}
void h_28(double *state, double *unused, double *out_1987303207216808992) {
   out_1987303207216808992[0] = state[5];
   out_1987303207216808992[1] = state[6];
}
void H_28(double *state, double *unused, double *out_961736427189882611) {
   out_961736427189882611[0] = 0;
   out_961736427189882611[1] = 0;
   out_961736427189882611[2] = 0;
   out_961736427189882611[3] = 0;
   out_961736427189882611[4] = 0;
   out_961736427189882611[5] = 1;
   out_961736427189882611[6] = 0;
   out_961736427189882611[7] = 0;
   out_961736427189882611[8] = 0;
   out_961736427189882611[9] = 0;
   out_961736427189882611[10] = 0;
   out_961736427189882611[11] = 0;
   out_961736427189882611[12] = 0;
   out_961736427189882611[13] = 0;
   out_961736427189882611[14] = 1;
   out_961736427189882611[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
