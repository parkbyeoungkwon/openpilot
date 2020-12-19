/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1451439341321623872);
void inv_err_fun(double *nom_x, double *true_x, double *out_344070343895654797);
void H_mod_fun(double *state, double *out_6601920486145661280);
void f_fun(double *state, double dt, double *out_2794436004131785345);
void F_fun(double *state, double dt, double *out_5440589488485531549);
void h_3(double *state, double *unused, double *out_8415310835275276442);
void H_3(double *state, double *unused, double *out_6009972283942761356);
void h_4(double *state, double *unused, double *out_5266440984635510722);
void H_4(double *state, double *unused, double *out_6355391271005691219);
void h_9(double *state, double *unused, double *out_5366006574141341390);
void H_9(double *state, double *unused, double *out_3013143561500021576);
void h_10(double *state, double *unused, double *out_4734333965803075872);
void H_10(double *state, double *unused, double *out_7735835397222591117);
void h_12(double *state, double *unused, double *out_7399726381137438920);
void H_12(double *state, double *unused, double *out_4248835759266481650);
void h_31(double *state, double *unused, double *out_6245193899023448090);
void H_31(double *state, double *unused, double *out_6923925534409901471);
void h_32(double *state, double *unused, double *out_2365153731292143974);
void H_32(double *state, double *unused, double *out_2317795955160842075);
void h_13(double *state, double *unused, double *out_3524022978074027593);
void H_13(double *state, double *unused, double *out_3872572480824100519);
void h_14(double *state, double *unused, double *out_5366006574141341390);
void H_14(double *state, double *unused, double *out_3013143561500021576);
void h_19(double *state, double *unused, double *out_14825494111016659);
void H_19(double *state, double *unused, double *out_7445823412602452757);
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