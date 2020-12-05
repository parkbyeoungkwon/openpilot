/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8479723950175235536);
void inv_err_fun(double *nom_x, double *true_x, double *out_7506973548091101796);
void H_mod_fun(double *state, double *out_3736768173062447715);
void f_fun(double *state, double dt, double *out_1303557177327667900);
void F_fun(double *state, double dt, double *out_8706792500591335821);
void h_25(double *state, double *unused, double *out_4714553408162648138);
void H_25(double *state, double *unused, double *out_954762370280144593);
void h_24(double *state, double *unused, double *out_2450204843206926825);
void H_24(double *state, double *unused, double *out_5134946454097058836);
void h_30(double *state, double *unused, double *out_4113135250549284584);
void H_30(double *state, double *unused, double *out_7878082792480223743);
void h_26(double *state, double *unused, double *out_4907408161891191858);
void H_26(double *state, double *unused, double *out_2296198503975792145);
void h_27(double *state, double *unused, double *out_8730517397886057833);
void H_27(double *state, double *unused, double *out_6590500804643598431);
void h_29(double *state, double *unused, double *out_8887704066237186978);
void H_29(double *state, double *unused, double *out_8452641406335759918);
void h_28(double *state, double *unused, double *out_2532851924621361964);
void H_28(double *state, double *unused, double *out_6293965141375997914);
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
