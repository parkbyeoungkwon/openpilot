/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_688280526805196633);
void inv_err_fun(double *nom_x, double *true_x, double *out_2805558480149713234);
void H_mod_fun(double *state, double *out_1867856875021459668);
void f_fun(double *state, double dt, double *out_4137902222323051559);
void F_fun(double *state, double dt, double *out_8216550828788730968);
void h_3(double *state, double *unused, double *out_183962545612576017);
void H_3(double *state, double *unused, double *out_8703255726107869327);
void h_4(double *state, double *unused, double *out_2734359870586169679);
void H_4(double *state, double *unused, double *out_8276716634326261519);
void h_9(double *state, double *unused, double *out_140674699463488053);
void H_9(double *state, double *unused, double *out_5954211090455405549);
void h_10(double *state, double *unused, double *out_8216205269947479839);
void H_10(double *state, double *unused, double *out_4519651043412688457);
void h_12(double *state, double *unused, double *out_3587706104945266578);
void H_12(double *state, double *unused, double *out_6170161122587051950);
void h_31(double *state, double *unused, double *out_3899283765689625939);
void H_31(double *state, double *unused, double *out_8845250897730471771);
void h_32(double *state, double *unused, double *out_2305182974926560217);
void H_32(double *state, double *unused, double *out_653443471437503972);
void h_13(double *state, double *unused, double *out_3653743223012753910);
void H_13(double *state, double *unused, double *out_1523046281694359924);
void h_14(double *state, double *unused, double *out_140674699463488053);
void H_14(double *state, double *unused, double *out_5954211090455405549);
void h_19(double *state, double *unused, double *out_7847435509852868596);
void H_19(double *state, double *unused, double *out_1521531239352974368);
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