#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "MatC.h"

// USER LEVEL
enum SOLVE_METHOD {
	EKF,
	SLSQ,
};
enum USING_OBS {
	CODE_ONLY,
	CODE_PHASE,
	PHASE_ONLY,
};
const SOLVE_METHOD METHOD = SLSQ;
USING_OBS USING = CODE_PHASE;
#define REF_PRN 28
#define WRONG_AMB_PRN 17
#define PAUSE_EPOCH 2241


// PROGRAM LEVEL
#define OBS_FILE "RemoteL1L2.obs"  
#define SAT_FILE "Satellites.sat"
#define STA_FILE "BaseL1L2.obs"

// RECIEVER LEVEL
#define GPS_SAT_NUM 32
#define CHANNEL_NUM 12

const double station_xyz[3]{ -1633489.379677252, -3651627.182503625, 4952481.619549180 };
const double base_xyz[3]{ -1625352.136567826, -3653483.734198109, 4953733.892847151 };

// OBSERVATION EVENTS
#define SAT_EVENT   unsigned char
#define SAT_LOW        0b00000001
#define SAT_OUT        0b00000010
#define SAT_IN         0b00000100
#define CODE_OUTLIER   0b00001000
#define CYCLE_SLIP_1   0b00010000
#define HIGH_DOP       0b00100000
#define LACK_SAT       0b01000000
#define CYCLE_SLIP_2   0b10000000

// CONDITION PARAMETERS
#define ELEV_THRES 0.10
#define ELEV_LOW   0.10
#define CYCSLIP_THRES_1 1 //cy
#define CYCSLIP_THRES_2 0.2 // m
#define GDOP_THRES 5.0
#define OUTLIER_THRES 20 // meters

// ESTIMATOR PARAMETERS
#define CODE_SIGMA0 8
#define PHASE_SIGMA0 0.08
#define FIXED_SIGMA0 100

// SLSQ
#define SLSQ_SYS_NOI 0.000000
#define FIXED_SYS_NOI 0.0
// EKF
#define EKF_SYS_NOI 0.1
#define FIRST_OBS_NOI 6e6
#define EKF_INTERVAL 1

// SOS PARAMETERS
#define SOS_THRES 3
#define SEARCH_SPACE_CONST 1.5
#define SOS_VALID_RES 0.50
#define SOS_START 1400

// GPS LEVEL
#define LAMBDA_1 0.190293672798365
#define LAMBDA_2 0.244210213424568

struct node {
	int epoch;
	int sat_num;
	SAT_EVENT ev;
	node * next;
};

struct perfect_segment {
	int start_time;
	int last_for;
	int sn;
};

struct obs_epoch {
	double PRN;
	double Time;
	double C1;
	double L1;
	double D1;
	double L2;
};

struct sat_epoch {
	double PRN;
	double Time;
	double P[3];
	double V[3];
};

struct gnss_ekf
{
	Matrix * X;
	Matrix * Xp;
	Matrix * Dx;
	Matrix * Dp;
	Matrix * F;

	Matrix * Ft;
	Matrix * Q;

	Matrix * Z;
	Matrix * H;
	Matrix * Dz;
	Matrix * Zp;

	Matrix * K;
	Matrix * V;
	Matrix * I;
};

void find_min(double * samples, int length, int & index, double & value, int except = -1)
{
	int start_index = 0;
	while (start_index == except) start_index++;

	index = start_index;
	value = samples[start_index];

	for (int i = start_index + 1; i < length; i++)
		if (samples[i] < value && i != except) {
			index = i; value = samples[i];
		}
}

double distance(const double * p1, const double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}

double elev(const double * sat_pos, const double * user_pos)
{

	double dpos[3] = { 0 };
	double ori[3]{ 0,0,0 };
	dpos[0] = sat_pos[0] - user_pos[0];
	dpos[1] = sat_pos[1] - user_pos[1];
	dpos[2] = sat_pos[2] - user_pos[2];

	double user_distance_to_earth = distance(user_pos, ori);

	double mod = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
	if (fabs(user_distance_to_earth * mod < 1.0)) {
		return M_PI_2;
	}
	else {
		double m = dpos[0] * user_pos[0] + dpos[1] * user_pos[1] + dpos[2] * user_pos[2];
		double n = m / (mod * user_distance_to_earth);
		return M_PI_2 - acos(n);
	}
}


obs_epoch OBS[CHANNEL_NUM];
sat_epoch SAT[CHANNEL_NUM];
obs_epoch STA[CHANNEL_NUM];

obs_epoch P_OBS[CHANNEL_NUM];
sat_epoch P_SAT[CHANNEL_NUM];
obs_epoch P_STA[CHANNEL_NUM];

FILE * ofp;
FILE * nfp;
FILE * sfp;

obs_epoch * S_OBS[CHANNEL_NUM];
sat_epoch * S_SAT[CHANNEL_NUM];
obs_epoch * S_STA[CHANNEL_NUM];

obs_epoch * P_S_OBS[CHANNEL_NUM];
sat_epoch * P_S_SAT[CHANNEL_NUM];
obs_epoch * P_S_STA[CHANNEL_NUM];

bool solve_available[GPS_SAT_NUM];
bool last_solve_available[GPS_SAT_NUM];

int sat_num;
int last_sat_num;

int current_time;

int history_amount = 0;
node * history = new node();
node * current = history;
SAT_EVENT ev;

// SOLUTION LEVEL
double sat_elev[GPS_SAT_NUM];
perfect_segment sigment;
double solution[3] = { -1625352.17084393, -3653483.75114927, 4953733.86925805 };
//double solution[3] = { -1633489.41308729 , -3651627.19449714, 4952481.59981920 };

double ambiguity[GPS_SAT_NUM] = {0};
// DOUBLE-DIFFER
int ref_sat = 0;
double ref_ele = 0;
Matrix * D = NULL;
int prn_table[GPS_SAT_NUM];
// SLSQ
Matrix * Q = NULL;
Matrix * Qinv = NULL;
// EKF
gnss_ekf ekf;

// SOS LEVEL
const int ref_amb[CHANNEL_NUM] = { 1, -12, -1, -10, -9, -20, -13, 8, -2 };
double sos_res_mag = 0;
double sigma_hat = 0;
int search_diameter[CHANNEL_NUM] = { 0 };
int int_n_essemble[CHANNEL_NUM + 1] = { 0 };
int ** possible_n;
int sos_search_num = 0;
bool sos_open = false;
Matrix ** int_n = NULL;
double * sos_omega = NULL;
int min_index = 0;
double sos_ratio = 0;
bool sos_fixed = false;
Matrix * cx = NULL;
Matrix * cn = NULL;
Matrix * n  = NULL;

// ANALYZE LEVEL
double phase_rate_value_1[GPS_SAT_NUM];
double phase_rate_value_2[GPS_SAT_NUM];
double phase_rate_value_3[GPS_SAT_NUM];
double phase_rate_value_4[GPS_SAT_NUM];
double phase_example;
double gdop;
double pseudorange_outliers[GPS_SAT_NUM];
double residuals[GPS_SAT_NUM];


void append_history(SAT_EVENT e, int sn)
{
	node * next = new node();
	next->epoch = current_time;
	next->ev = e;
	next->sat_num = sn;
	next->next = NULL;
	current->next = next;
	current = current->next;
	history_amount++;
}

// phase rate method.
inline double cycle_slip_1(obs_epoch * current, obs_epoch * last)
{
	return fabs(last->L1 + (current->D1 + last->D1) / 2 - current->L1);
}

// phase difference method.
inline double cycle_slip_2(obs_epoch * current, obs_epoch * last)
{
	return fabs((current->L1 - last->L1)* LAMBDA_1 - (current->L2 - last->L2) * LAMBDA_2);
}

void ekf_create(double DeltaT, int sat_num)
{
	double ttd2 = 0.5 * DeltaT * DeltaT;
	double tttd6 = DeltaT * DeltaT * DeltaT / 6;
	int num_state = 2 + sat_num;
	ekf.X = malloc_mat(num_state, 1);
	ekf.X->data[0][0] = station_xyz[0];
	ekf.X->data[1][0] = station_xyz[1];
	ekf.X->data[2][0] = station_xyz[2];


	ekf.Xp = malloc_mat(num_state, 1);
	ekf.Dx = eyes(num_state, FIRST_OBS_NOI);
	ekf.Dp = eyes(num_state);
	ekf.I = eyes(num_state);

	ekf.Z = NULL;
	ekf.Dz = NULL;
	ekf.H = NULL;
	ekf.Ft = NULL;
	ekf.K = NULL;
	ekf.V = NULL;
	ekf.Zp = NULL;

	ekf.F = eyes(num_state);
	ekf.Q = eyes(num_state, EKF_SYS_NOI);

	mat_trans(ekf.F, ekf.Ft);
}

void ekf_predict()
{
	mat_multiply(ekf.F, ekf.X, ekf.Xp);
	Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL;
	mat_multiply(ekf.F, ekf.Dx, temp1);
	mat_multiply(temp1, ekf.Ft, temp2);
	mat_sum(temp2, ekf.Q, ekf.Dp);

	free_mat(temp1);
	free_mat(temp2);
	free_mat(temp3);
	free_mat(temp4);
}

void ekf_fetch_observation()
{
	free_mat(ekf.H);
	free_mat(ekf.Z);
	free_mat(ekf.Dz);
	free_mat(ekf.Zp);

	ekf.H = malloc_mat(sat_num * 2, 2 + sat_num);
	ekf.Z = malloc_mat(sat_num * 2, 1);
	ekf.Dz = malloc_mat(sat_num * 2, sat_num * 2);
	ekf.Zp = malloc_mat(sat_num * 2, 1);

	double x[3] = {
		ekf.Xp->data[0][0],
		ekf.Xp->data[1][0],
		ekf.Xp->data[2][0]
	};

	double * DX0 = (double*)alloca(sat_num * sizeof(double));
	double * DY0 = (double*)alloca(sat_num * sizeof(double));
	double * DZ0 = (double*)alloca(sat_num * sizeof(double));
	double *   S = (double*)alloca(sat_num * sizeof(double));
	double *  S2 = (double*)alloca(sat_num * sizeof(double));

	for (int i = 0; i < sat_num; i++)
	{
		DX0[i] = S_SAT[i]->P[0] - x[0];
		DY0[i] = S_SAT[i]->P[1] - x[1];
		DZ0[i] = S_SAT[i]->P[2] - x[2];
		S[i] = sqrt(DX0[i] * DX0[i] + DY0[i] * DY0[i] + DZ0[i] * DZ0[i]);

		S2[i] = distance(base_xyz, S_SAT[i]->P);
	}
	for (int j = 0; j < sat_num; j++)
	{
		int i_c = j * 2;     // code
		int i_p = j * 2 + 1; // phase

		double sinE = sin(sat_elev[(int)(S_SAT[j]->PRN - 1)]);
		ekf.Dz->data[i_c][i_c] = CODE_SIGMA0 * CODE_SIGMA0 / sinE / sinE;
		ekf.Dz->data[i_p][i_p] = PHASE_SIGMA0 * PHASE_SIGMA0 / sinE / sinE;

		ekf.Z->data[i_c][0] = S_OBS[j]->C1 - S_STA[j]->C1 + S2[j];
		ekf.Z->data[i_p][0] = S_OBS[j]->L1 * LAMBDA_1 - S_STA[j]->L1 * LAMBDA_1 + S2[j];

		ekf.H->data[i_c][0] = -DX0[j] / S[j];
		ekf.H->data[i_c][1] = -DY0[j] / S[j];
		ekf.H->data[i_c][2] = -DZ0[j] / S[j];

		ekf.H->data[i_p][0] = -DX0[j] / S[j];
		ekf.H->data[i_p][1] = -DY0[j] / S[j];
		ekf.H->data[i_p][2] = -DZ0[j] / S[j];

		ekf.Zp->data[i_c][0] = S[j];
		ekf.Zp->data[i_p][0] = S[j];
	}

	Matrix * La = ekf.Z;
	Matrix * Lpa = ekf.Zp;
	Matrix * Aa = ekf.H;
	Matrix * Cla = ekf.Dz, *temp = NULL, *Dt = NULL;

	ekf.H = NULL;
	ekf.Z = NULL;
	ekf.Zp = NULL;
	ekf.Dz = NULL;

	mat_multiply(D, La, ekf.Z);
	mat_multiply(D, Lpa, ekf.Zp);
	mat_multiply(D, Aa, ekf.H);

	mat_multiply(D, Cla, temp);
	mat_trans(D, Dt);
	mat_multiply(temp, Dt, ekf.Dz);

	free_mat(Aa);
	free_mat(La);
	free_mat(Cla);
	free_mat(temp);
	free_mat(Dt);

	for (int i = 0; i < sat_num - 1; i++)
	{
		int i_p = i * 2 + 1; // phase
		ekf.H->data[i_p][i + 3] = LAMBDA_1;
		ekf.Zp->data[i_p][0] += ekf.Xp->data[i + 3][0] * LAMBDA_1;
	}
}

void ekf_execute()
{
	//Kalman Gain
	Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL;
	Matrix * Ht = NULL;
	mat_trans(ekf.H, Ht);
	mat_multiply(ekf.Dp, Ht, temp1);
	mat_multiply(ekf.H, ekf.Dp, temp2);
	mat_multiply(temp2, Ht, temp3);
	mat_sum(temp3, ekf.Dz);
	mat_inv(temp3, temp4);
	mat_multiply(temp1, temp4, ekf.K);
	free_mat(temp1); free_mat(temp2); free_mat(temp3); free_mat(temp4);

	//Kalman Innovation
	mat_minus(ekf.Z, ekf.Zp, ekf.V);

	//State Update
	mat_multiply(ekf.K, ekf.V, temp1);
	mat_sum(ekf.Xp, temp1, ekf.X);

	//Covariance Matrix
	mat_multiply(ekf.K, ekf.H, temp2);
	mat_minus(ekf.I, temp2, temp3);
	mat_multiply(temp3, ekf.Dp, ekf.Dx);

	free_mat(temp1); free_mat(temp2); free_mat(temp3);
	free_mat(ekf.K);
}

bool ekf_solve()
{
	ekf_predict();
	ekf_fetch_observation();
	ekf_execute();

	solution[0] = ekf.X->data[0][0];
	solution[1] = ekf.X->data[1][0];
	solution[2] = ekf.X->data[2][0];
	for (int i = 0; i < sat_num - 1; i++)
		ambiguity[prn_table[i] - 1] = ekf.X->data[i + 3][0];
	
	return true;
}

void inverse_phase(obs_epoch * ptr)
{
	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		ptr[i].L1 = -ptr[i].L1;
		ptr[i].L2 = -ptr[i].L2;
	}
}

bool outlier_check()
{
	for (int i = 0; i < sat_num; i++)
	{
		double dis_ba = distance(S_SAT[i]->P, base_xyz);
		double dis_st = distance(S_SAT[i]->P, station_xyz);
		double obs = S_OBS[i]->C1 - S_STA[i]->C1;
		double act = dis_st - dis_ba;
		pseudorange_outliers[(int)(S_OBS[i]->PRN - 1)] = fabs(obs - act);
		if (fabs(obs - act) >= OUTLIER_THRES) 
			return false;
	}
	return true;
}

bool dop_check()
{
	double * DX0 = (double*)alloca(sat_num * sizeof(double));
	double * DY0 = (double*)alloca(sat_num * sizeof(double));
	double * DZ0 = (double*)alloca(sat_num * sizeof(double)); 
	double *   S = (double*)alloca(sat_num * sizeof(double));  

	Matrix * A = malloc_mat(sat_num, 4);  

	for (int j = 0; j < sat_num; j++)
	{
		DX0[j] = S_SAT[j]->P[0] - station_xyz[0];
		DY0[j] = S_SAT[j]->P[1] - station_xyz[1];
		DZ0[j] = S_SAT[j]->P[2] - station_xyz[2];
		S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);
	}
	for (int j = 0; j < sat_num; j++)
	{
		A->data[j][0] = -DX0[j] / S[j]; 
		A->data[j][1] = -DY0[j] / S[j]; 
		A->data[j][2] = -DZ0[j] / S[j]; 
		A->data[j][3] = 1;               
	}

	Matrix * At = NULL, *AtA = NULL, *Q = NULL;
	mat_trans(A, At);
	mat_multiply(At, A, AtA);
	mat_inv(AtA, Q);

	gdop = sqrt(
		Q->data[0][0] + Q->data[1][1] + Q->data[2][2] + Q->data[3][3]);

	free_mat(A); free_mat(At); free_mat(Q); free_mat(AtA);
	if (gdop >= GDOP_THRES)
		return false;
	return gdop < GDOP_THRES;
}

void overall_reset()
{
	last_sat_num = sat_num;
	sat_num = 0;
	ev = 0;
	gdop = 0;

	memcpy(last_solve_available, solve_available, sizeof(bool) * GPS_SAT_NUM);
	memset(solve_available, 0, sizeof(bool) * GPS_SAT_NUM);

	memset(phase_rate_value_1, 0, sizeof(double) * GPS_SAT_NUM);
	memset(phase_rate_value_2, 0, sizeof(double) * GPS_SAT_NUM);
	memset(phase_rate_value_3, 0, sizeof(double) * GPS_SAT_NUM);
	memset(phase_rate_value_4, 0, sizeof(double) * GPS_SAT_NUM);
	memset(residuals, 0, sizeof(double) * GPS_SAT_NUM);
	memset(pseudorange_outliers, 0, sizeof(double) * GPS_SAT_NUM);

	memset(sat_elev, 0, sizeof(double) * GPS_SAT_NUM);

	//memset(prn_table, 0, sizeof(int) * GPS_SAT_NUM);

	memcpy(P_OBS, OBS, sizeof(obs_epoch) * CHANNEL_NUM);
	memcpy(P_SAT, SAT, sizeof(sat_epoch) * CHANNEL_NUM);
	memcpy(P_STA, STA, sizeof(obs_epoch) * CHANNEL_NUM);
}

bool overall_check()
{
	current_time = (int)OBS->Time;
	

	for (int i = 0; i < CHANNEL_NUM; i++)
	{
		int prn = (int)OBS[i].PRN;
		if (prn == 0)continue;
		for (int j = 0; j < CHANNEL_NUM; j++)
		{
			if (SAT[j].PRN == prn)
			{
				for (int k = 0; k < CHANNEL_NUM; k++)
				{
					if (STA[k].PRN == prn)
					{
						for (int l = 0; l < CHANNEL_NUM; l++)
						{
							if (P_OBS[l].PRN == prn)
							{
								for (int m = 0; m < CHANNEL_NUM; m++)
								{
									if (P_SAT[m].PRN == prn)
									{
										for (int n = 0; n < CHANNEL_NUM; n++)
										{
											if (P_STA[n].PRN == prn)
											{
												// elevation
												sat_elev[prn - 1] = elev(SAT[j].P, station_xyz);
												if (sat_elev[prn - 1] <= ELEV_THRES)
													goto next;
												if( sat_elev[prn - 1] <= ELEV_LOW)
													ev |= SAT_LOW;
												
												// cycleslip
												phase_rate_value_1[prn - 1] = cycle_slip_1(OBS + i, P_OBS + l);
												phase_rate_value_2[prn - 1] = cycle_slip_2(OBS + i, P_OBS + l);
												phase_rate_value_3[prn - 1] = cycle_slip_1(STA + k, P_STA + n);
												phase_rate_value_4[prn - 1] = cycle_slip_2(STA + k, P_STA + n);

												if (phase_rate_value_1[prn - 1] >= CYCSLIP_THRES_1)
												{
													ev |= CYCLE_SLIP_1;
												}
												if (phase_rate_value_2[prn - 1] >= CYCSLIP_THRES_2)
												{
													ev |= CYCLE_SLIP_2;
												}

												if (phase_rate_value_3[prn - 1] >= CYCSLIP_THRES_1)
												{
													ev |= CYCLE_SLIP_1;
												}
												if (phase_rate_value_4[prn - 1] >= CYCSLIP_THRES_2)
												{
													ev |= CYCLE_SLIP_2;
												}

												if (!last_solve_available[prn - 1])
												{
													ev |= SAT_IN;
													last_solve_available[prn - 1] = true;
												}

												solve_available[prn - 1] = true;
												S_OBS[sat_num] = OBS + i;
												S_SAT[sat_num] = SAT + j;
												S_STA[sat_num] = STA + k;
												P_S_OBS[sat_num] = P_OBS + l;
												P_S_SAT[sat_num] = P_SAT + m;
												P_S_STA[sat_num] = P_STA + n;
												sat_num++;
												goto next;
											}
										}
										goto next;
									}
								}
								goto next;
							}
						}
						goto next;
					}
				}
				goto next;
			}
		}
	next:;
	}

	if (memcmp(last_solve_available, solve_available, sizeof(bool) * GPS_SAT_NUM) != 0)
		ev |= SAT_OUT;
	if (!outlier_check())
		ev |= CODE_OUTLIER;
	if (sat_num < 5)
		return !(ev |= LACK_SAT);
	if (!dop_check())
		ev |= HIGH_DOP;

	return true;
}

bool find_perfect_sigment()
{
	for (node * i = history->next; i != NULL; i = i->next)
	{
	next:;
		if (i->ev != 0) continue;
		for (node * j = i; j != NULL; j = j->next)
		{
			if (j->ev != 0)
			{
				if (sigment.last_for < j->epoch - 1 - i->epoch) {
					sigment.start_time = i->epoch;
					sigment.last_for = j->epoch - 1 - i->epoch;
					sigment.sn = i->sat_num;
				}
				i = j->next;
				
				goto next;
			}
		}
	}

	return true;
}

void fetch_ref_sat(int prn = 0)
{
	if (prn <= 0) {
		ref_sat = 0;
		ref_ele = sat_elev[0];

		for (int i = 1; i < sat_num; i++)
		{
			if (ref_ele < sat_elev[(int)(S_SAT[i]->PRN) - 1])
			{
				ref_ele = sat_elev[(int)(S_SAT[i]->PRN) - 1];
				ref_sat = i;
			}
		}
	}
	else
	{
		for (int i = 0; i < sat_num; i++)
		{
			if (S_SAT[i]->PRN == prn)
			{
				ref_sat = i;
				ref_ele = sat_elev[prn - 1];
				goto next;
			}
		}

	}
	throw - 1;
next:;

	for (int i = 0; i < sat_num - 1; i++)
	{
		prn_table[i] = (int)S_SAT[i + (i >= ref_sat ? 1 : 0)]->PRN;
	}
	int row = 0, col = 0;
	switch (USING) {
	case CODE_ONLY:
	case PHASE_ONLY:
		row = sat_num - 1;
		col = sat_num;
		D = malloc_mat(row, col);
		for (int i = 0; i < row; i++)
		{
			D->data[i][i + (i >= ref_sat ? 1 : 0)] = 1;
			D->data[i][ref_sat] = -1;
		}
		break;
	case CODE_PHASE:
		row = 2 * sat_num - 2;
		col = 2 * sat_num;
		D = malloc_mat(row, col);
		for (int i = 0; i < row; i += 2)
		{
			D->data[i][i + (i >= ref_sat * 2 ? 2 : 0)] = 1;
			D->data[i + 1][i + (i >= ref_sat * 2 ? 3 : 1)] = 1;
			D->data[i][ref_sat * 2] = -1;
			D->data[i + 1][ref_sat * 2 + 1] = -1;
		}
		break;
	default:
		throw "under construction";
	}
}

bool slsq_solve()
{

	Matrix * L = NULL;
	Matrix * A = NULL;
	Matrix * Cl = NULL;
	Matrix * r = NULL;
	Matrix * δ = NULL;
	Matrix * Cx = NULL;

	switch (USING)
	{
	case CODE_ONLY:
	case PHASE_ONLY:
		L = malloc_mat(sat_num, 1);
		A = malloc_mat(sat_num, 3);
		Cl = malloc_mat(sat_num, sat_num);
		r = malloc_mat(sat_num, 1);
		δ = malloc_mat(3, 1);
		Cx = malloc_mat(3, 3);
		break;
	case CODE_PHASE:
		L = malloc_mat(sat_num * 2, 1);
		A = malloc_mat(sat_num * 2, sat_num + 2);
		Cl = malloc_mat(sat_num * 2, sat_num * 2);
		r = malloc_mat(sat_num * 2, 1);
		δ = malloc_mat(sat_num + 2, 1);
		Cx = malloc_mat(sat_num + 2, sat_num + 2);
		break;
	default:
		throw "under construction";
	}

	double * DX0 = (double*)alloca(sat_num * sizeof(double));
	double * DY0 = (double*)alloca(sat_num * sizeof(double));
	double * DZ0 = (double*)alloca(sat_num * sizeof(double));
	double * S = (double*)alloca(sat_num * sizeof(double)); 
	double * S2 = (double*)alloca(sat_num * sizeof(double)); 


	//for (int i = 0; i < LS_MAX_ITER; i++) 
	{
		for (int j = 0; j < sat_num; j++)
		{
			DX0[j] = S_SAT[j]->P[0] - solution[0];
			DY0[j] = S_SAT[j]->P[1] - solution[1];
			DZ0[j] = S_SAT[j]->P[2] - solution[2];
			S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);

			S2[j] = distance(base_xyz, S_SAT[j]->P);
		}

		if(USING == CODE_PHASE)
			for (int j = 0; j < sat_num; j++)
			{
				int i_c = j * 2;     // code
				int i_p = j * 2 + 1; // phase

				double sinE = sin(sat_elev[(int)(S_SAT[j]->PRN - 1)]);
				Cl->data[i_c][i_c] = CODE_SIGMA0 * CODE_SIGMA0 / sinE / sinE;
				Cl->data[i_p][i_p] = PHASE_SIGMA0 * PHASE_SIGMA0 / sinE / sinE;

				L->data[i_c][0] = S_OBS[j]->C1 - S_STA[j]->C1 + S2[j] - S[j];
				L->data[i_p][0] = (S_OBS[j]->L1 - S_STA[j]->L1) * LAMBDA_1 + S2[j] - S[j];

				A->data[i_c][0] = -DX0[j] / S[j];
				A->data[i_c][1] = -DY0[j] / S[j];
				A->data[i_c][2] = -DZ0[j] / S[j];

				A->data[i_p][0] = -DX0[j] / S[j];
				A->data[i_p][1] = -DY0[j] / S[j];
				A->data[i_p][2] = -DZ0[j] / S[j];
			}
		else if(USING == CODE_ONLY)
			for (int j = 0; j < sat_num; j++)
			{
				double sinE = sin(sat_elev[(int)(S_SAT[j]->PRN - 1)]);
				Cl->data[j][j] = CODE_SIGMA0 * CODE_SIGMA0 / sinE / sinE;
				L->data[j][0] = S_OBS[j]->C1 - S_STA[j]->C1 + S2[j] - S[j];
				A->data[j][0] = -DX0[j] / S[j];
				A->data[j][1] = -DY0[j] / S[j];
				A->data[j][2] = -DZ0[j] / S[j];
			}
		else if (USING == PHASE_ONLY)
			for (int j = 0; j < sat_num; j++)
			{
				double sinE = sin(sat_elev[(int)(S_SAT[j]->PRN - 1)]);
				Cl->data[j][j] = FIXED_SIGMA0 * PHASE_SIGMA0 * PHASE_SIGMA0 / sinE / sinE;
				L->data[j][0] = (S_OBS[j]->L1 - S_STA[j]->L1) * LAMBDA_1 + S2[j] - S[j];
				A->data[j][0] = -DX0[j] / S[j];
				A->data[j][1] = -DY0[j] / S[j];
				A->data[j][2] = -DZ0[j] / S[j];
			}


		Matrix * La = L;
		Matrix * Aa = A;
		Matrix * Cla = Cl, *temp = NULL, *Dt = NULL;

		A = NULL;
		L = NULL;
		Cl = NULL;
		mat_multiply(D, La, L);
		mat_multiply(D, Aa, A); 

		mat_multiply(D, Cla, temp);
		mat_trans(D, Dt);
		mat_multiply(temp, Dt, Cl);

		free_mat(Aa);
		free_mat(La);
		free_mat(Cla);
		free_mat(temp);
		free_mat(Dt);

		if (USING == CODE_PHASE)
			for (int i = 0; i < sat_num - 1; i++)
			{
				int i_p = i * 2 + 1; // phase
				A->data[i_p][i + 3] = LAMBDA_1;
				L->data[i_p][0] -= ambiguity[prn_table[i] - 1] * LAMBDA_1;
			}
		else if (USING == PHASE_ONLY)
			for (int i = 0; i < sat_num - 1; i++)
				L->data[i][0] -= ambiguity[prn_table[i] - 1] * LAMBDA_1;

		SLS(L, A, Cl, δ, Qinv, Q);

		Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL, *temp5 = NULL;
		Matrix * temp6 = NULL;
		Matrix * temp7 = NULL, *temp8 = NULL, *temp9 = NULL, *temp10 = NULL, * temp11 = NULL;
		Matrix * temp12 = NULL, *temp13 = NULL, *cn_inv = NULL;
		free_mat(cx); free_mat(cn);

		mat_multiply(A, δ, temp1);
		mat_minus(L, temp1, temp2);
		mat_trans(temp2, temp3);
		mat_inv(Cl, temp4);
		mat_multiply(temp3, temp4, temp5);
		mat_multiply(temp5, temp2, temp6);
		sos_res_mag = temp6->data[0][0];

		sigma_hat = /*USING == PHASE_ONLY ? sos_res_mag / (sat_num - 4) :*/ 1;//sos_res_mag / (sat_num - 4);
		mat_multiply(Q, sigma_hat, cx);

		if (USING == PHASE_ONLY)
		{
			for (int i = 0; i < sat_num - 1; i++)
			{
				residuals[prn_table[i] - 1] = temp2->data[i][0];
			}
		}

		//sub_mat(Q, 3, sat_num + 1, 3, sat_num + 1, temp7);
		//sub_mat(temp2, 3, sat_num + 1, 0, 0, temp8);
		//mat_trans(temp8, temp9);
		//mat_multiply(temp9, temp4, temp10);
		//mat_multiply(temp10, temp4, temp11);
		//mat_multiply(temp7, temp11->data[0][0] / sat_num - 7, cn);
		
		//mat_output(cx, "cx");


		
		//cn = malloc_mat(sat_num - 1, sat_num - 1);
		//for (int i = 0; i < sat_num - 1; i++)
		//	cn->data[i][i] = cx->data[i + 3][i + 3];

		if (USING == CODE_PHASE) {
			sub_mat(cx, 3, sat_num + 1, 3, sat_num + 1, cn);
			for (int j = 0; j < sat_num - 1; j++)
			{
				ambiguity[prn_table[j] - 1] += δ->data[j + 3][0];
				Q->data[j + 3][j + 3] += SLSQ_SYS_NOI;
			}
		}

		if (sos_open && USING == CODE_PHASE)
		{
			Matrix * cn2 = NULL;
			mat_multiply(cn, 1, cn2);
			mat_inv(cn2, cn_inv);

			for (int i = 0; i < sat_num - 1; i++)
				n->data[i][0] = ambiguity[prn_table[i] - 1];


			for (int i = 0; i < sos_search_num; i++)
			{
				Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL, *temp5 = NULL;
				
				mat_minus(n, int_n[i], temp1);
				mat_trans(temp1, temp2);
				
				mat_multiply(temp2, cn_inv, temp4);
				mat_multiply(temp4, temp1, temp5);
				/*mat_multiply(temp2, temp1, temp5);*/
				sos_omega[i] += temp5->data[0][0];
				free_mat(temp1);
				free_mat(temp2);
				free_mat(temp3);
				free_mat(temp4);
				free_mat(temp5);
			}
		}
		for (int j = 0; j < 3; j++) {
			solution[j] = solution[j] + δ->data[j][0];
			Q->data[j][j] += USING == PHASE_ONLY ? FIXED_SYS_NOI : SLSQ_SYS_NOI;
		}
		mat_inv(Q, Qinv);
		free_mat(temp1); free_mat(temp2); free_mat(temp3);
		free_mat(temp4); free_mat(temp5); free_mat(temp6);
		free_mat(temp7); free_mat(temp8); free_mat(temp9);
		free_mat(temp10); free_mat(temp11); free_mat(temp12);
		free_mat(temp13); free_mat(cn_inv);
		

	}
	free_mat(L); free_mat(A); free_mat(Cl);
	free_mat(r); free_mat(δ);
	free_mat(Cx);
	return true;
}

void n_enum(int sat_index)
{
	if (sat_index == sat_num) {
		int_n[sos_search_num] = malloc_mat(sat_num - 1, 1);
		for (int i = 0; i < sat_num - 1; i++)
			int_n[sos_search_num]->data[i][0] = int_n_essemble[i];
		//mat_output(int_n[sos_search_num], "int_n");
		sos_search_num++;
		return;
	}

	int possibility = sat_index == sat_num - 1 ? 1 : search_diameter[sat_index];
	for (int i = 0; i < possibility; i++)
	{
		int_n_essemble[sat_index] = possible_n[sat_index][i];
		n_enum(sat_index + 1);
	}
}

void start_sos()
{
	printf("==========START SOS PROCESS==========\n");
	sos_open = true;
	sos_search_num = 1;
	possible_n = new int*[sat_num];
	possible_n[sat_num - 1] = new int[CHANNEL_NUM];
	int mark = 9;
	for (int i = 0; i < sat_num - 1; i++)
	{
		bool has = false;

		double n = ambiguity[prn_table[i] - 1];
		double r = sqrt(cn->data[i][i]) * SEARCH_SPACE_CONST;
		int upper = (int)round(n + r);
		int litter = (int)round(n - r);
		possible_n[i] = new int[upper - litter + 1];
		for (int j = litter; j <= upper; j++) {
			possible_n[i][j - litter] = j;
			if (ref_amb[i] == j) has = true;
			printf("%d,", j);
		}
		search_diameter[i] = upper - litter + 1;
		sos_search_num *= (search_diameter[i]);
		printf("%s\n", has ? "yes" : "no");
		if (has)mark--;
	}
	printf("remain: %d, space: %d\n", mark, sos_search_num);
	int_n = new Matrix*[sos_search_num];
	n = malloc_mat(sat_num - 1, 1);
	for (int i = 0; i < sat_num - 1; i++)
		n->data[i][0] = ambiguity[prn_table[i] - 1];
	sos_omega = new double[sos_search_num];
	memset(sos_omega, 0, sizeof(double) * sos_search_num);

	sos_search_num = 0;
	printf("=====================================\n");
	n_enum(0);
}

void check_sos()
{
	double min = 0;
	double min2 = 0;
	
	int min2_index = 0;

	//FILE * fp = fopen("omega.txt", "w");

	//fclose(fp);

	find_min(sos_omega, sos_search_num, min_index, min);
	find_min(sos_omega, sos_search_num, min2_index, min2, min_index);
	
	sos_ratio = min2 / min;
	if (!sos_fixed && sos_ratio > SOS_THRES) {
		sos_fixed = true;
		printf("======= SOS ALGORITHM FIXED  ========\n");
		printf("=====================================\n");
	}
}

FILE * ccsl;
FILE * ccsl2;
FILE * ccsl3;
FILE * ccsl4;
FILE * elevs;
FILE * outl;
FILE * gdopf;
FILE * eve;
FILE * cxf;
FILE * omf;

FILE * exam;
FILE * outfp;
FILE * outfp_fw;
FILE * outfp_f;
FILE * cxff;
FILE * resif;
FILE * resifw;
FILE * amb;

void output_1()
{
	for (int i = 0; i < GPS_SAT_NUM; i++)
	{
		fprintf(ccsl, "%lf\t", phase_rate_value_1[i]);
		fprintf(ccsl2, "%lf\t", phase_rate_value_2[i]);
		fprintf(ccsl3, "%lf\t", phase_rate_value_3[i]);
		fprintf(ccsl4, "%lf\t", phase_rate_value_4[i]);
		fprintf(elevs, "%lf\t", sat_elev[i]);
		fprintf(outl, "%lf\t", pseudorange_outliers[i]);
	}
	fprintf(ccsl, "\n");
	fprintf(ccsl2, "\n");
	fprintf(ccsl3, "\n");
	fprintf(ccsl4, "\n");
	fprintf(elevs, "\n");
	fprintf(outl, "\n");

	fprintf(gdopf, "%lf\n", gdop);
	
	static unsigned char and[8]{
		0b00000001,
		0b00000010,
		0b00000100,
		0b00001000,
		0b00010000,
		0b00100000,
		0b01000000,
		0b10000000
	};
	for (int i = 0; i < 8; i++)
	{
		if (ev != 0)
		{
			int i = 1;
		}
		fprintf(eve, "%d\t", (((ev & and[i]) >> (i)) == 1) ? 1 : 0);
	}
	fprintf(eve, "\n");

	fprintf(exam, "%lf\n", OBS[0].L1);
}

void output_2()
{
	fprintf(outfp, "%lf\t%lf\t%lf\n", solution[0], solution[1], solution[2]);

	for (int i = 0; i < GPS_SAT_NUM; i++) {
		fprintf(amb, "%lf\t", ambiguity[i]);
	}

	for (int i = 0; i < sos_search_num; i++)
	{
		fprintf(omf, "%lf\t", sos_omega[i]);
	}
	fprintf(omf, "\n");

	for (int i = 0; i < cx->rows; i++)
	{
		fprintf(cxf, "%lf\t", cx->data[i][i]);
	}
	fprintf(cxf, "\n");

	fprintf(amb, "\n");
}

void output_3()
{
	fprintf(outfp_f, "%lf\t%lf\t%lf\n", solution[0], solution[1], solution[2]);
	for (int i = 0; i < cx->rows; i++)
	{
		fprintf(cxff, "%lf\t", cx->data[i][i]);
	}
	for (int i = 0; i < GPS_SAT_NUM; i++) {
		fprintf(resif, "%lf\t", residuals[i]);
	}
	fprintf(resif, "\n");
	fprintf(cxff, "\n");
}

void output_4()
{
	for (int i = 0; i < GPS_SAT_NUM; i++) {
		fprintf(resifw, "%lf\t", residuals[i]);
	}
	fprintf(resifw, "\n");
	fprintf(outfp_fw, "%lf\t%lf\t%lf\n", solution[0], solution[1], solution[2]);
}


void check_observation()
{
	ccsl = fopen("ccsl.txt", "w");
	ccsl2 = fopen("ccsl2.txt", "w");
	ccsl3 = fopen("ccsl3.txt", "w");
	ccsl4 = fopen("ccsl4.txt", "w");
	elevs = fopen("elevs.txt", "w");
	outl = fopen("outl.txt", "w");
	exam = fopen("exam.txt", "w");
	gdopf = fopen("gdopf.txt", "w");
	eve = fopen("eve.txt", "w");

	ofp = fopen(OBS_FILE, "rb");
	nfp = fopen(SAT_FILE, "rb");
	sfp = fopen(STA_FILE, "rb");
	int epoch_num = 1;
	while (!feof(ofp) && !feof(nfp) && !feof(sfp))
	{
		if (!fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp)) break;
		if (!fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp)) break;
		if (!fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp)) break;
		inverse_phase(OBS);
		inverse_phase(STA);
		if (epoch_num == PAUSE_EPOCH)
		{
			int j = 0;
		}
		overall_check();
		append_history(ev, sat_num);

		output_1();

		overall_reset();
		epoch_num++;
	}

	_fcloseall();
}

bool fix_n()
{
	sat_num = int_n[0]->rows + 1;
	for (int i = 0; i < sat_num - 1; i++)
	{
		if (round(ambiguity[prn_table[i] - 1]) != int_n[min_index]->data[i][0] && int_n[min_index]->data[i][0] != ref_amb[i])
			return false;
		
		ambiguity[prn_table[i] - 1] = int_n[min_index]->data[i][0];
	}
	return true;
}

void solve_fixed_wrong()
{
	USING = PHASE_ONLY;
	outfp_fw = fopen("outfw.txt", "w");
	resifw = fopen("resiw.txt", "w");
	if (!fix_n())throw - 1;

	// make it wrong
	ambiguity[WRONG_AMB_PRN - 1] --;

	ofp = fopen(OBS_FILE, "rb");
	nfp = fopen(SAT_FILE, "rb");
	sfp = fopen(STA_FILE, "rb");
	

	if (METHOD == EKF)
		ekf_create(EKF_INTERVAL, sigment.sn);
	else if (METHOD == SLSQ)
	{
		free_mat(Q);
		free_mat(Qinv);
	}

	while (!feof(ofp) && !feof(nfp) && !feof(sfp))
	{
		fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp);
		fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp);
		fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp);

		inverse_phase(OBS);
		inverse_phase(STA);
		overall_check();

		if (current_time - sigment.start_time == PAUSE_EPOCH)
		{
			int j = 0;
		}
		if (OBS->Time >= sigment.start_time && OBS->Time < sigment.last_for + sigment.start_time)
		{
			fetch_ref_sat(REF_PRN);
			if (METHOD == SLSQ)
				slsq_solve();
			else if (METHOD == EKF)
				ekf_solve();
			output_4();
		}
		overall_reset();
	}

	_fcloseall();
}

void solve_fixed()
{
	USING = PHASE_ONLY;
	outfp_f = fopen("outf.txt", "w");
	cxff = fopen("cxff.txt", "w");
	resif = fopen("resi.txt", "w");
	
	if (!fix_n())throw - 1;

	ofp = fopen(OBS_FILE, "rb");
	nfp = fopen(SAT_FILE, "rb");
	sfp = fopen(STA_FILE, "rb");

	if (METHOD == EKF)
		ekf_create(EKF_INTERVAL, sigment.sn);
	else if (METHOD == SLSQ)
	{
		free_mat(Q);
		free_mat(Qinv);
	}

	while (!feof(ofp) && !feof(nfp) && !feof(sfp))
	{
		fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp);
		fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp);
		fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp);

		inverse_phase(OBS);
		inverse_phase(STA);
		overall_check();

		if (current_time - sigment.start_time == PAUSE_EPOCH)
		{
			int j = 0;
		}
		if (OBS->Time >= sigment.start_time && OBS->Time < sigment.last_for + sigment.start_time)
		{
			fetch_ref_sat(REF_PRN);
			if (METHOD == SLSQ)
				slsq_solve();
			else if (METHOD == EKF)
				ekf_solve();
			output_3();
		}
		overall_reset();
	}

	_fcloseall();
}

void solve_float()
{
	outfp = fopen("out.txt", "w");
	ofp = fopen(OBS_FILE, "rb");
	nfp = fopen(SAT_FILE, "rb");
	sfp = fopen(STA_FILE, "rb");
	omf = fopen("ome.txt", "w");
	amb = fopen("amb.txt", "w");
	cxf = fopen("cxf.txt", "w");
	if (METHOD == EKF)
		ekf_create(EKF_INTERVAL, sigment.sn);
	else if (METHOD == SLSQ)
	{
		free_mat(Q);
		free_mat(Qinv);
	}

	while (!feof(ofp) && !feof(nfp) && !feof(sfp))
	{
		

		fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp);
		fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp);
		fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp);

		inverse_phase(OBS);
		inverse_phase(STA);
		overall_check();

		if (current_time - sigment.start_time == PAUSE_EPOCH)
		{
			int j = 0;
		}
		if (OBS->Time >= sigment.start_time && OBS->Time < sigment.last_for + sigment.start_time)
		{
			
			fetch_ref_sat(REF_PRN);

			if (METHOD == SLSQ)
				slsq_solve();
			else if (METHOD == EKF)
				ekf_solve();
			if(sos_open)
				check_sos();
			output_2();
		}

		if (current_time - sigment.start_time == SOS_START)
			start_sos();

		overall_reset();
	}

	_fcloseall();
}

int main()
{
	if (USING == CODE_ONLY && METHOD == EKF) return 0;

	printf("=====================================\n");
	printf("====== SEARCHING OBSERVATION  =======\n");
	check_observation();
	printf("=====================================\n");
	printf("========= LOCATING SIGMENT  =========\n");
	
	find_perfect_sigment();
	printf("=====================================\n");
	printf("=========== SOLVING FLOAT ===========\n");
	solve_float();

	if (sos_fixed)
	{
		printf("=========== SOLVING FIXED ===========\n");
		solve_fixed();
		printf("=====================================\n");
		printf("=========== SOLVING WRONG ===========\n");
		solve_fixed_wrong();
	}
	printf("=====================================\n");
	system("pause");
}
