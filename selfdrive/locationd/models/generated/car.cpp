
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
void err_fun(double *nom_x, double *delta_x, double *out_7198846295097160569) {
   out_7198846295097160569[0] = delta_x[0] + nom_x[0];
   out_7198846295097160569[1] = delta_x[1] + nom_x[1];
   out_7198846295097160569[2] = delta_x[2] + nom_x[2];
   out_7198846295097160569[3] = delta_x[3] + nom_x[3];
   out_7198846295097160569[4] = delta_x[4] + nom_x[4];
   out_7198846295097160569[5] = delta_x[5] + nom_x[5];
   out_7198846295097160569[6] = delta_x[6] + nom_x[6];
   out_7198846295097160569[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_11262323420528530) {
   out_11262323420528530[0] = -nom_x[0] + true_x[0];
   out_11262323420528530[1] = -nom_x[1] + true_x[1];
   out_11262323420528530[2] = -nom_x[2] + true_x[2];
   out_11262323420528530[3] = -nom_x[3] + true_x[3];
   out_11262323420528530[4] = -nom_x[4] + true_x[4];
   out_11262323420528530[5] = -nom_x[5] + true_x[5];
   out_11262323420528530[6] = -nom_x[6] + true_x[6];
   out_11262323420528530[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1318688685977315366) {
   out_1318688685977315366[0] = 1.0;
   out_1318688685977315366[1] = 0.0;
   out_1318688685977315366[2] = 0.0;
   out_1318688685977315366[3] = 0.0;
   out_1318688685977315366[4] = 0.0;
   out_1318688685977315366[5] = 0.0;
   out_1318688685977315366[6] = 0.0;
   out_1318688685977315366[7] = 0.0;
   out_1318688685977315366[8] = 0.0;
   out_1318688685977315366[9] = 1.0;
   out_1318688685977315366[10] = 0.0;
   out_1318688685977315366[11] = 0.0;
   out_1318688685977315366[12] = 0.0;
   out_1318688685977315366[13] = 0.0;
   out_1318688685977315366[14] = 0.0;
   out_1318688685977315366[15] = 0.0;
   out_1318688685977315366[16] = 0.0;
   out_1318688685977315366[17] = 0.0;
   out_1318688685977315366[18] = 1.0;
   out_1318688685977315366[19] = 0.0;
   out_1318688685977315366[20] = 0.0;
   out_1318688685977315366[21] = 0.0;
   out_1318688685977315366[22] = 0.0;
   out_1318688685977315366[23] = 0.0;
   out_1318688685977315366[24] = 0.0;
   out_1318688685977315366[25] = 0.0;
   out_1318688685977315366[26] = 0.0;
   out_1318688685977315366[27] = 1.0;
   out_1318688685977315366[28] = 0.0;
   out_1318688685977315366[29] = 0.0;
   out_1318688685977315366[30] = 0.0;
   out_1318688685977315366[31] = 0.0;
   out_1318688685977315366[32] = 0.0;
   out_1318688685977315366[33] = 0.0;
   out_1318688685977315366[34] = 0.0;
   out_1318688685977315366[35] = 0.0;
   out_1318688685977315366[36] = 1.0;
   out_1318688685977315366[37] = 0.0;
   out_1318688685977315366[38] = 0.0;
   out_1318688685977315366[39] = 0.0;
   out_1318688685977315366[40] = 0.0;
   out_1318688685977315366[41] = 0.0;
   out_1318688685977315366[42] = 0.0;
   out_1318688685977315366[43] = 0.0;
   out_1318688685977315366[44] = 0.0;
   out_1318688685977315366[45] = 1.0;
   out_1318688685977315366[46] = 0.0;
   out_1318688685977315366[47] = 0.0;
   out_1318688685977315366[48] = 0.0;
   out_1318688685977315366[49] = 0.0;
   out_1318688685977315366[50] = 0.0;
   out_1318688685977315366[51] = 0.0;
   out_1318688685977315366[52] = 0.0;
   out_1318688685977315366[53] = 0.0;
   out_1318688685977315366[54] = 1.0;
   out_1318688685977315366[55] = 0.0;
   out_1318688685977315366[56] = 0.0;
   out_1318688685977315366[57] = 0.0;
   out_1318688685977315366[58] = 0.0;
   out_1318688685977315366[59] = 0.0;
   out_1318688685977315366[60] = 0.0;
   out_1318688685977315366[61] = 0.0;
   out_1318688685977315366[62] = 0.0;
   out_1318688685977315366[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9057978945416650967) {
   out_9057978945416650967[0] = state[0];
   out_9057978945416650967[1] = state[1];
   out_9057978945416650967[2] = state[2];
   out_9057978945416650967[3] = state[3];
   out_9057978945416650967[4] = state[4];
   out_9057978945416650967[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9057978945416650967[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9057978945416650967[7] = state[7];
}
void F_fun(double *state, double dt, double *out_5819831454079858889) {
   out_5819831454079858889[0] = 1;
   out_5819831454079858889[1] = 0;
   out_5819831454079858889[2] = 0;
   out_5819831454079858889[3] = 0;
   out_5819831454079858889[4] = 0;
   out_5819831454079858889[5] = 0;
   out_5819831454079858889[6] = 0;
   out_5819831454079858889[7] = 0;
   out_5819831454079858889[8] = 0;
   out_5819831454079858889[9] = 1;
   out_5819831454079858889[10] = 0;
   out_5819831454079858889[11] = 0;
   out_5819831454079858889[12] = 0;
   out_5819831454079858889[13] = 0;
   out_5819831454079858889[14] = 0;
   out_5819831454079858889[15] = 0;
   out_5819831454079858889[16] = 0;
   out_5819831454079858889[17] = 0;
   out_5819831454079858889[18] = 1;
   out_5819831454079858889[19] = 0;
   out_5819831454079858889[20] = 0;
   out_5819831454079858889[21] = 0;
   out_5819831454079858889[22] = 0;
   out_5819831454079858889[23] = 0;
   out_5819831454079858889[24] = 0;
   out_5819831454079858889[25] = 0;
   out_5819831454079858889[26] = 0;
   out_5819831454079858889[27] = 1;
   out_5819831454079858889[28] = 0;
   out_5819831454079858889[29] = 0;
   out_5819831454079858889[30] = 0;
   out_5819831454079858889[31] = 0;
   out_5819831454079858889[32] = 0;
   out_5819831454079858889[33] = 0;
   out_5819831454079858889[34] = 0;
   out_5819831454079858889[35] = 0;
   out_5819831454079858889[36] = 1;
   out_5819831454079858889[37] = 0;
   out_5819831454079858889[38] = 0;
   out_5819831454079858889[39] = 0;
   out_5819831454079858889[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5819831454079858889[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5819831454079858889[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5819831454079858889[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5819831454079858889[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5819831454079858889[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5819831454079858889[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5819831454079858889[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5819831454079858889[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5819831454079858889[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5819831454079858889[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5819831454079858889[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5819831454079858889[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5819831454079858889[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5819831454079858889[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5819831454079858889[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5819831454079858889[56] = 0;
   out_5819831454079858889[57] = 0;
   out_5819831454079858889[58] = 0;
   out_5819831454079858889[59] = 0;
   out_5819831454079858889[60] = 0;
   out_5819831454079858889[61] = 0;
   out_5819831454079858889[62] = 0;
   out_5819831454079858889[63] = 1;
}
void h_25(double *state, double *unused, double *out_6721412496810265817) {
   out_6721412496810265817[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3084851460078232124) {
   out_3084851460078232124[0] = 0;
   out_3084851460078232124[1] = 0;
   out_3084851460078232124[2] = 0;
   out_3084851460078232124[3] = 0;
   out_3084851460078232124[4] = 0;
   out_3084851460078232124[5] = 0;
   out_3084851460078232124[6] = 1;
   out_3084851460078232124[7] = 0;
}
void h_24(double *state, double *unused, double *out_6691078089342205072) {
   out_6691078089342205072[0] = state[4];
   out_6691078089342205072[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4347240524273972256) {
   out_4347240524273972256[0] = 0;
   out_4347240524273972256[1] = 0;
   out_4347240524273972256[2] = 0;
   out_4347240524273972256[3] = 0;
   out_4347240524273972256[4] = 1;
   out_4347240524273972256[5] = 0;
   out_4347240524273972256[6] = 0;
   out_4347240524273972256[7] = 0;
   out_4347240524273972256[8] = 0;
   out_4347240524273972256[9] = 0;
   out_4347240524273972256[10] = 0;
   out_4347240524273972256[11] = 0;
   out_4347240524273972256[12] = 0;
   out_4347240524273972256[13] = 1;
   out_4347240524273972256[14] = 0;
   out_4347240524273972256[15] = 0;
}
void h_30(double *state, double *unused, double *out_1041852760010406624) {
   out_1041852760010406624[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6529047450870951156) {
   out_6529047450870951156[0] = 0;
   out_6529047450870951156[1] = 0;
   out_6529047450870951156[2] = 0;
   out_6529047450870951156[3] = 0;
   out_6529047450870951156[4] = 1;
   out_6529047450870951156[5] = 0;
   out_6529047450870951156[6] = 0;
   out_6529047450870951156[7] = 0;
}
void h_26(double *state, double *unused, double *out_4447374810821714131) {
   out_4447374810821714131[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8789444615017441397) {
   out_8789444615017441397[0] = 0;
   out_8789444615017441397[1] = 0;
   out_8789444615017441397[2] = 0;
   out_8789444615017441397[3] = 0;
   out_8789444615017441397[4] = 0;
   out_8789444615017441397[5] = 0;
   out_8789444615017441397[6] = 0;
   out_8789444615017441397[7] = 1;
}
void h_27(double *state, double *unused, double *out_6275742451873912680) {
   out_6275742451873912680[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7816629438707576468) {
   out_7816629438707576468[0] = 0;
   out_7816629438707576468[1] = 0;
   out_7816629438707576468[2] = 0;
   out_7816629438707576468[3] = 1;
   out_7816629438707576468[4] = 0;
   out_7816629438707576468[5] = 0;
   out_7816629438707576468[6] = 0;
   out_7816629438707576468[7] = 0;
}
void h_29(double *state, double *unused, double *out_2861974832809938502) {
   out_2861974832809938502[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8093897853709768507) {
   out_8093897853709768507[0] = 0;
   out_8093897853709768507[1] = 1;
   out_8093897853709768507[2] = 0;
   out_8093897853709768507[3] = 0;
   out_8093897853709768507[4] = 0;
   out_8093897853709768507[5] = 0;
   out_8093897853709768507[6] = 0;
   out_8093897853709768507[7] = 0;
}
void h_28(double *state, double *unused, double *out_2133854971428204598) {
   out_2133854971428204598[0] = state[5];
   out_2133854971428204598[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7081671071199084494) {
   out_7081671071199084494[0] = 0;
   out_7081671071199084494[1] = 0;
   out_7081671071199084494[2] = 0;
   out_7081671071199084494[3] = 0;
   out_7081671071199084494[4] = 0;
   out_7081671071199084494[5] = 1;
   out_7081671071199084494[6] = 0;
   out_7081671071199084494[7] = 0;
   out_7081671071199084494[8] = 0;
   out_7081671071199084494[9] = 0;
   out_7081671071199084494[10] = 0;
   out_7081671071199084494[11] = 0;
   out_7081671071199084494[12] = 0;
   out_7081671071199084494[13] = 0;
   out_7081671071199084494[14] = 1;
   out_7081671071199084494[15] = 0;
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
