/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7198846295097160569);
void inv_err_fun(double *nom_x, double *true_x, double *out_11262323420528530);
void H_mod_fun(double *state, double *out_1318688685977315366);
void f_fun(double *state, double dt, double *out_9057978945416650967);
void F_fun(double *state, double dt, double *out_5819831454079858889);
void h_25(double *state, double *unused, double *out_6721412496810265817);
void H_25(double *state, double *unused, double *out_3084851460078232124);
void h_24(double *state, double *unused, double *out_6691078089342205072);
void H_24(double *state, double *unused, double *out_4347240524273972256);
void h_30(double *state, double *unused, double *out_1041852760010406624);
void H_30(double *state, double *unused, double *out_6529047450870951156);
void h_26(double *state, double *unused, double *out_4447374810821714131);
void H_26(double *state, double *unused, double *out_8789444615017441397);
void h_27(double *state, double *unused, double *out_6275742451873912680);
void H_27(double *state, double *unused, double *out_7816629438707576468);
void h_29(double *state, double *unused, double *out_2861974832809938502);
void H_29(double *state, double *unused, double *out_8093897853709768507);
void h_28(double *state, double *unused, double *out_2133854971428204598);
void H_28(double *state, double *unused, double *out_7081671071199084494);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
