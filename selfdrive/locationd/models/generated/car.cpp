
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
void err_fun(double *nom_x, double *delta_x, double *out_8479723950175235536) {
   out_8479723950175235536[0] = delta_x[0] + nom_x[0];
   out_8479723950175235536[1] = delta_x[1] + nom_x[1];
   out_8479723950175235536[2] = delta_x[2] + nom_x[2];
   out_8479723950175235536[3] = delta_x[3] + nom_x[3];
   out_8479723950175235536[4] = delta_x[4] + nom_x[4];
   out_8479723950175235536[5] = delta_x[5] + nom_x[5];
   out_8479723950175235536[6] = delta_x[6] + nom_x[6];
   out_8479723950175235536[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7506973548091101796) {
   out_7506973548091101796[0] = -nom_x[0] + true_x[0];
   out_7506973548091101796[1] = -nom_x[1] + true_x[1];
   out_7506973548091101796[2] = -nom_x[2] + true_x[2];
   out_7506973548091101796[3] = -nom_x[3] + true_x[3];
   out_7506973548091101796[4] = -nom_x[4] + true_x[4];
   out_7506973548091101796[5] = -nom_x[5] + true_x[5];
   out_7506973548091101796[6] = -nom_x[6] + true_x[6];
   out_7506973548091101796[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3736768173062447715) {
   out_3736768173062447715[0] = 1.0;
   out_3736768173062447715[1] = 0.0;
   out_3736768173062447715[2] = 0.0;
   out_3736768173062447715[3] = 0.0;
   out_3736768173062447715[4] = 0.0;
   out_3736768173062447715[5] = 0.0;
   out_3736768173062447715[6] = 0.0;
   out_3736768173062447715[7] = 0.0;
   out_3736768173062447715[8] = 0.0;
   out_3736768173062447715[9] = 1.0;
   out_3736768173062447715[10] = 0.0;
   out_3736768173062447715[11] = 0.0;
   out_3736768173062447715[12] = 0.0;
   out_3736768173062447715[13] = 0.0;
   out_3736768173062447715[14] = 0.0;
   out_3736768173062447715[15] = 0.0;
   out_3736768173062447715[16] = 0.0;
   out_3736768173062447715[17] = 0.0;
   out_3736768173062447715[18] = 1.0;
   out_3736768173062447715[19] = 0.0;
   out_3736768173062447715[20] = 0.0;
   out_3736768173062447715[21] = 0.0;
   out_3736768173062447715[22] = 0.0;
   out_3736768173062447715[23] = 0.0;
   out_3736768173062447715[24] = 0.0;
   out_3736768173062447715[25] = 0.0;
   out_3736768173062447715[26] = 0.0;
   out_3736768173062447715[27] = 1.0;
   out_3736768173062447715[28] = 0.0;
   out_3736768173062447715[29] = 0.0;
   out_3736768173062447715[30] = 0.0;
   out_3736768173062447715[31] = 0.0;
   out_3736768173062447715[32] = 0.0;
   out_3736768173062447715[33] = 0.0;
   out_3736768173062447715[34] = 0.0;
   out_3736768173062447715[35] = 0.0;
   out_3736768173062447715[36] = 1.0;
   out_3736768173062447715[37] = 0.0;
   out_3736768173062447715[38] = 0.0;
   out_3736768173062447715[39] = 0.0;
   out_3736768173062447715[40] = 0.0;
   out_3736768173062447715[41] = 0.0;
   out_3736768173062447715[42] = 0.0;
   out_3736768173062447715[43] = 0.0;
   out_3736768173062447715[44] = 0.0;
   out_3736768173062447715[45] = 1.0;
   out_3736768173062447715[46] = 0.0;
   out_3736768173062447715[47] = 0.0;
   out_3736768173062447715[48] = 0.0;
   out_3736768173062447715[49] = 0.0;
   out_3736768173062447715[50] = 0.0;
   out_3736768173062447715[51] = 0.0;
   out_3736768173062447715[52] = 0.0;
   out_3736768173062447715[53] = 0.0;
   out_3736768173062447715[54] = 1.0;
   out_3736768173062447715[55] = 0.0;
   out_3736768173062447715[56] = 0.0;
   out_3736768173062447715[57] = 0.0;
   out_3736768173062447715[58] = 0.0;
   out_3736768173062447715[59] = 0.0;
   out_3736768173062447715[60] = 0.0;
   out_3736768173062447715[61] = 0.0;
   out_3736768173062447715[62] = 0.0;
   out_3736768173062447715[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_1303557177327667900) {
   out_1303557177327667900[0] = state[0];
   out_1303557177327667900[1] = state[1];
   out_1303557177327667900[2] = state[2];
   out_1303557177327667900[3] = state[3];
   out_1303557177327667900[4] = state[4];
   out_1303557177327667900[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1303557177327667900[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1303557177327667900[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8706792500591335821) {
   out_8706792500591335821[0] = 1;
   out_8706792500591335821[1] = 0;
   out_8706792500591335821[2] = 0;
   out_8706792500591335821[3] = 0;
   out_8706792500591335821[4] = 0;
   out_8706792500591335821[5] = 0;
   out_8706792500591335821[6] = 0;
   out_8706792500591335821[7] = 0;
   out_8706792500591335821[8] = 0;
   out_8706792500591335821[9] = 1;
   out_8706792500591335821[10] = 0;
   out_8706792500591335821[11] = 0;
   out_8706792500591335821[12] = 0;
   out_8706792500591335821[13] = 0;
   out_8706792500591335821[14] = 0;
   out_8706792500591335821[15] = 0;
   out_8706792500591335821[16] = 0;
   out_8706792500591335821[17] = 0;
   out_8706792500591335821[18] = 1;
   out_8706792500591335821[19] = 0;
   out_8706792500591335821[20] = 0;
   out_8706792500591335821[21] = 0;
   out_8706792500591335821[22] = 0;
   out_8706792500591335821[23] = 0;
   out_8706792500591335821[24] = 0;
   out_8706792500591335821[25] = 0;
   out_8706792500591335821[26] = 0;
   out_8706792500591335821[27] = 1;
   out_8706792500591335821[28] = 0;
   out_8706792500591335821[29] = 0;
   out_8706792500591335821[30] = 0;
   out_8706792500591335821[31] = 0;
   out_8706792500591335821[32] = 0;
   out_8706792500591335821[33] = 0;
   out_8706792500591335821[34] = 0;
   out_8706792500591335821[35] = 0;
   out_8706792500591335821[36] = 1;
   out_8706792500591335821[37] = 0;
   out_8706792500591335821[38] = 0;
   out_8706792500591335821[39] = 0;
   out_8706792500591335821[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8706792500591335821[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8706792500591335821[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8706792500591335821[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8706792500591335821[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8706792500591335821[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8706792500591335821[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8706792500591335821[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8706792500591335821[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8706792500591335821[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8706792500591335821[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8706792500591335821[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8706792500591335821[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8706792500591335821[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8706792500591335821[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8706792500591335821[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8706792500591335821[56] = 0;
   out_8706792500591335821[57] = 0;
   out_8706792500591335821[58] = 0;
   out_8706792500591335821[59] = 0;
   out_8706792500591335821[60] = 0;
   out_8706792500591335821[61] = 0;
   out_8706792500591335821[62] = 0;
   out_8706792500591335821[63] = 1;
}
void h_25(double *state, double *unused, double *out_4714553408162648138) {
   out_4714553408162648138[0] = state[6];
}
void H_25(double *state, double *unused, double *out_954762370280144593) {
   out_954762370280144593[0] = 0;
   out_954762370280144593[1] = 0;
   out_954762370280144593[2] = 0;
   out_954762370280144593[3] = 0;
   out_954762370280144593[4] = 0;
   out_954762370280144593[5] = 0;
   out_954762370280144593[6] = 1;
   out_954762370280144593[7] = 0;
}
void h_24(double *state, double *unused, double *out_2450204843206926825) {
   out_2450204843206926825[0] = state[4];
   out_2450204843206926825[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5134946454097058836) {
   out_5134946454097058836[0] = 0;
   out_5134946454097058836[1] = 0;
   out_5134946454097058836[2] = 0;
   out_5134946454097058836[3] = 0;
   out_5134946454097058836[4] = 1;
   out_5134946454097058836[5] = 0;
   out_5134946454097058836[6] = 0;
   out_5134946454097058836[7] = 0;
   out_5134946454097058836[8] = 0;
   out_5134946454097058836[9] = 0;
   out_5134946454097058836[10] = 0;
   out_5134946454097058836[11] = 0;
   out_5134946454097058836[12] = 0;
   out_5134946454097058836[13] = 1;
   out_5134946454097058836[14] = 0;
   out_5134946454097058836[15] = 0;
}
void h_30(double *state, double *unused, double *out_4113135250549284584) {
   out_4113135250549284584[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7878082792480223743) {
   out_7878082792480223743[0] = 0;
   out_7878082792480223743[1] = 0;
   out_7878082792480223743[2] = 0;
   out_7878082792480223743[3] = 0;
   out_7878082792480223743[4] = 1;
   out_7878082792480223743[5] = 0;
   out_7878082792480223743[6] = 0;
   out_7878082792480223743[7] = 0;
}
void h_26(double *state, double *unused, double *out_4907408161891191858) {
   out_4907408161891191858[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2296198503975792145) {
   out_2296198503975792145[0] = 0;
   out_2296198503975792145[1] = 0;
   out_2296198503975792145[2] = 0;
   out_2296198503975792145[3] = 0;
   out_2296198503975792145[4] = 0;
   out_2296198503975792145[5] = 0;
   out_2296198503975792145[6] = 0;
   out_2296198503975792145[7] = 1;
}
void h_27(double *state, double *unused, double *out_8730517397886057833) {
   out_8730517397886057833[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6590500804643598431) {
   out_6590500804643598431[0] = 0;
   out_6590500804643598431[1] = 0;
   out_6590500804643598431[2] = 0;
   out_6590500804643598431[3] = 1;
   out_6590500804643598431[4] = 0;
   out_6590500804643598431[5] = 0;
   out_6590500804643598431[6] = 0;
   out_6590500804643598431[7] = 0;
}
void h_29(double *state, double *unused, double *out_8887704066237186978) {
   out_8887704066237186978[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8452641406335759918) {
   out_8452641406335759918[0] = 0;
   out_8452641406335759918[1] = 1;
   out_8452641406335759918[2] = 0;
   out_8452641406335759918[3] = 0;
   out_8452641406335759918[4] = 0;
   out_8452641406335759918[5] = 0;
   out_8452641406335759918[6] = 0;
   out_8452641406335759918[7] = 0;
}
void h_28(double *state, double *unused, double *out_2532851924621361964) {
   out_2532851924621361964[0] = state[5];
   out_2532851924621361964[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6293965141375997914) {
   out_6293965141375997914[0] = 0;
   out_6293965141375997914[1] = 0;
   out_6293965141375997914[2] = 0;
   out_6293965141375997914[3] = 0;
   out_6293965141375997914[4] = 0;
   out_6293965141375997914[5] = 1;
   out_6293965141375997914[6] = 0;
   out_6293965141375997914[7] = 0;
   out_6293965141375997914[8] = 0;
   out_6293965141375997914[9] = 0;
   out_6293965141375997914[10] = 0;
   out_6293965141375997914[11] = 0;
   out_6293965141375997914[12] = 0;
   out_6293965141375997914[13] = 0;
   out_6293965141375997914[14] = 1;
   out_6293965141375997914[15] = 0;
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
