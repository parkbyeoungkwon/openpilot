/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6676735503722724430);
void inv_err_fun(double *nom_x, double *true_x, double *out_3440255382924592903);
void H_mod_fun(double *state, double *out_6857355617626461041);
void f_fun(double *state, double dt, double *out_6643389448991676382);
void F_fun(double *state, double dt, double *out_8607567293629774541);
void h_3(double *state, double *unused, double *out_2559027154466927835);
void H_3(double *state, double *unused, double *out_1797205890141238734);
void h_4(double *state, double *unused, double *out_3410191088690539276);
void H_4(double *state, double *unused, double *out_7793698294824612578);
void h_9(double *state, double *unused, double *out_3986247643202311009);
void H_9(double *state, double *unused, double *out_5471192750953756608);
void h_10(double *state, double *unused, double *out_5781822902471576951);
void H_10(double *state, double *unused, double *out_2452826951470157496);
void h_12(double *state, double *unused, double *out_4836164021828696522);
void H_12(double *state, double *unused, double *out_5687142783085403009);
void h_31(double *state, double *unused, double *out_5799404479008061282);
void H_31(double *state, double *unused, double *out_8362232558228822830);
void h_32(double *state, double *unused, double *out_2900044907023652711);
void H_32(double *state, double *unused, double *out_8932023104681401314);
void h_13(double *state, double *unused, double *out_1442324608208524788);
void H_13(double *state, double *unused, double *out_4376924510876197122);
void h_14(double *state, double *unused, double *out_3986247643202311009);
void H_14(double *state, double *unused, double *out_5471192750953756608);
void h_19(double *state, double *unused, double *out_8714221944303198677);
void H_19(double *state, double *unused, double *out_1038512899851325427);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);