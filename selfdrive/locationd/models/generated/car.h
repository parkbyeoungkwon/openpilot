/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6439155458356320514);
void inv_err_fun(double *nom_x, double *true_x, double *out_3744964575875050616);
void H_mod_fun(double *state, double *out_6996129398019384578);
void f_fun(double *state, double dt, double *out_8382255420543571790);
void F_fun(double *state, double dt, double *out_3259366650936592814);
void h_25(double *state, double *unused, double *out_2795798451157236819);
void H_25(double *state, double *unused, double *out_3511382771814294907);
void h_24(double *state, double *unused, double *out_3100482985859878334);
void H_24(double *state, double *unused, double *out_7979568905426377477);
void h_30(double *state, double *unused, double *out_7887680365731642356);
void H_30(double *state, double *unused, double *out_5321462390946073429);
void h_26(double *state, double *unused, double *out_3123920273549208447);
void H_26(double *state, double *unused, double *out_2193210383124914366);
void h_27(double *state, double *unused, double *out_8394801905458863088);
void H_27(double *state, double *unused, double *out_4033880403109448117);
void h_29(double *state, double *unused, double *out_5557390440317063615);
void H_29(double *state, double *unused, double *out_1497663621817241476);
void h_28(double *state, double *unused, double *out_1987303207216808992);
void H_28(double *state, double *unused, double *out_961736427189882611);
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
