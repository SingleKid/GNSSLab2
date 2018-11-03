#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "MatC.h"
#include "cmn.h"

// USER LEVEL
enum SOLVE_METHOD {
	EKF,
	SLSQ,
};
const SOLVE_METHOD METHOD = EKF;

// PROGRAM LEVEL
#define OBS_FILE "RemoteL1L2.obs"  
#define SAT_FILE "Satellites.sat"
#define STA_FILE "BaseL1L2.obs"

// RECIEVER LEVEL
#define GPS_SAT_NUM 32
#define CHANNEL_NUM 12
const double station_xyz[3]{ -1633489.41308729 , -3651627.19449714, 4952481.59981920 };
const double base_xyz[3]{ -1625352.17084393, -3653483.75114927, 4953733.86925805 };

// OBSERVATION EVENTS
#define SAT_EVENT unsigned char
#define SAT_LOW        0b00000001
#define SAT_OUT        0b00000010
#define SAT_IN         0b00000100
#define CODE_OUTLIER   0b00001000
#define CYCLE_SLIP_1   0b00010000
#define HIGH_DOP       0b00100000
#define LACK_SAT       0b01000000
#define CYCLE_SLIP_2   0b10000000

// CONDITION PARAMETERS
#define ELEV_THRES 0.1
#define ELEV_LOW   0.1
#define CYCSLIP_THRES_1 1 //cy
#define CYCSLIP_THRES_2 0.2 // m
#define GDOP_THRES 5
#define OUTLIER_THRES 20 // meters

// ESTIMATOR PARAMETERS
#define OBS_SIGMA0 1
#define SLSQ_SYS_NOI 0.5
#define EKF_SYS_NOI 0.1
#define FIRST_OBS_NOI 6e6
#define EKF_INTERVAL 1

// GPS LEVEL
#define LAMBDA_1 0.190293672798365
#define LAMBDA_2 0.244210213424568

struct node {
	int epoch;
	SAT_EVENT ev;
	node * next;
};

struct perfect_segment {
	int start_time;
	int last_for;
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
	Matrix * T;

	Matrix * Ft;
	Matrix * Tt;
	Matrix * De;

	Matrix * Z;
	Matrix * H;
	Matrix * Dz;
	Matrix * Zp;

	Matrix * K;
	Matrix * V;
};

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
double solution[3] = { 0,0,0 };
// DOUBLE-DIFFER
int ref_sat = 0;
double ref_ele = 0;
Matrix * D = NULL;
// SLSQ
Matrix * Q = NULL;
// EKF
gnss_ekf ekf;

// ANALYZE LEVEL
double phase_rate_value_1[GPS_SAT_NUM];
double phase_rate_value_2[GPS_SAT_NUM];

void append_history(SAT_EVENT e)
{
	node * next = new node();
	next->epoch = current_time;
	next->ev = e;
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

void ekf_create(double DeltaT)
{
	double ttd2 = 0.5 * DeltaT * DeltaT;
	double tttd6 = DeltaT * DeltaT * DeltaT / 6;
	ekf.X = malloc_mat(9, 1);
	ekf.Xp = malloc_mat(9, 1);
	ekf.Dx = eyes(9, FIRST_OBS_NOI);
	ekf.Dp = eyes(9);

	ekf.Z = NULL;
	ekf.Dz = NULL;
	ekf.H = NULL;
	ekf.Ft = NULL;
	ekf.Tt = NULL;
	ekf.K = NULL;
	ekf.V = NULL;
	ekf.Zp = NULL;

	ekf.F = malloc_mat(9, 9);
	for (int i = 0; i < 3; i++)
	{
		int offset = i * 3;
		ekf.F->data[offset][offset] = 1;
		ekf.F->data[offset][offset + 1] = DeltaT;
		ekf.F->data[offset][offset + 2] = ttd2;
		ekf.F->data[offset + 1][offset + 1] = 1;
		ekf.F->data[offset + 1][offset + 2] = DeltaT;
		ekf.F->data[offset + 2][offset + 2] = 1;
	}

	ekf.T = malloc_mat(9, 3);
	for (int i = 0; i < 3; i++)
	{
		ekf.T->data[3 * i][i] = tttd6;
		ekf.T->data[3 * i + 1][i] = ttd2;
		ekf.T->data[3 * i + 2][i] = DeltaT;
	}

	ekf.De = eyes(3, EKF_SYS_NOI);

	mat_trans(ekf.F, ekf.Ft);
	mat_trans(ekf.T, ekf.Tt);
}

void ekf_predict()
{
	mat_multiply(ekf.F, ekf.X, ekf.Xp);
	Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL;
	mat_multiply(ekf.F, ekf.Dx, temp1);
	mat_multiply(temp1, ekf.Ft, temp2);
	mat_multiply(ekf.T, ekf.De, temp3);
	mat_multiply(temp3, ekf.Tt, temp4);
	mat_sum(temp2, temp4, ekf.Dp);

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

	ekf.H = malloc_mat(sat_num, 9);
	ekf.Z = malloc_mat(sat_num, 1);
	ekf.Dz = malloc_mat(sat_num, sat_num);
	ekf.Zp = malloc_mat(sat_num, 1);

	double x[3] = {
		ekf.Xp->data[0][0],
		ekf.Xp->data[3][0],
		ekf.Xp->data[6][0]
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
		double sinE = sin(sat_elev[(int)(S_SAT[j]->PRN - 1)]);
		ekf.Z->data[j][0] = S_OBS[j]->C1 - S_STA[j]->C1 + S2[j];
		ekf.Dz->data[j][j] = OBS_SIGMA0 * OBS_SIGMA0 / sinE / sinE;

		ekf.H->data[j][0] = -DX0[j] / S[j];
		ekf.H->data[j][3] = -DY0[j] / S[j];
		ekf.H->data[j][6] = -DZ0[j] / S[j];

		ekf.Zp->data[j][0] = S[j];
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
	Matrix * I = eyes(9);
	mat_minus(I, temp2, temp3);
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
	solution[1] = ekf.X->data[3][0];
	solution[2] = ekf.X->data[6][0];
	
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

	double GDOP = sqrt(
		Q->data[0][0] + Q->data[1][1] + Q->data[2][2] + Q->data[3][3]);

	free_mat(A); free_mat(At); free_mat(Q); free_mat(AtA);
	if (GDOP >= GDOP_THRES) 
		return false;
	return GDOP < GDOP_THRES;
}

void overall_reset()
{
	last_sat_num = sat_num;
	sat_num = 0;
	ev = 0;
	memcpy(last_solve_available, solve_available, sizeof(bool) * GPS_SAT_NUM);
	memset(solve_available, 0, sizeof(bool) * GPS_SAT_NUM);

	memset(phase_rate_value_1, 0, sizeof(double) * GPS_SAT_NUM);
	memset(phase_rate_value_2, 0, sizeof(double) * GPS_SAT_NUM);
	memset(sat_elev, 0, sizeof(double) * GPS_SAT_NUM);

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
												if (phase_rate_value_1[prn - 1] >= CYCSLIP_THRES_1)
												{
													ev |= CYCLE_SLIP_1;
												}
												if (phase_rate_value_2[prn - 1] >= CYCSLIP_THRES_2)
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
				}
				i = j->next;
				
				goto next;
			}
		}
	}

	return true;
}

void fetch_ref_sat()
{
	ref_sat = 0;
	ref_ele = sat_elev[0];

	for (int i = 1; i < sat_num; i++)
	{
		if (ref_ele < sat_elev[(int)(S_SAT[i]->PRN)])
		{
			ref_ele = sat_elev[(int)(S_SAT[i]->PRN)];
			ref_sat = i;
		}
	}

	D = malloc_mat(sat_num - 1, sat_num);
	for (int i = 0; i < sat_num - 1; i++)
	{
		D->data[i][i + (i >= ref_sat ? 1 : 0)] = 1;
		D->data[i][ref_sat] = -1;
	}
}

bool slsq_solve()
{
	Matrix * L = malloc_mat(sat_num, 1);
	Matrix * A = malloc_mat(sat_num, 3);
	Matrix * Cl = malloc_mat(sat_num, sat_num);
	Matrix * r = malloc_mat(sat_num, 1);
	Matrix * δ = malloc_mat(3, 1);
	Matrix * Cx = malloc_mat(3, 3);



	double * DX0 = (double*)alloca(sat_num * sizeof(double));
	double * DY0 = (double*)alloca(sat_num * sizeof(double));
	double * DZ0 = (double*)alloca(sat_num * sizeof(double));
	double * S = (double*)alloca(sat_num * sizeof(double)); 
	double * S2 = (double*)alloca(sat_num * sizeof(double)); 

	double last_solution[4] = { 0,0,0,0 };

	//for (int i = 0; i < LS_MAX_ITER; i++) 
	{
		memcpy(last_solution, solution, sizeof(double) * 4); 
		for (int j = 0; j < sat_num; j++)
		{
			DX0[j] = S_SAT[j]->P[0] - solution[0];
			DY0[j] = S_SAT[j]->P[1] - solution[1];
			DZ0[j] = S_SAT[j]->P[2] - solution[2];
			S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);

			S2[j] = distance(base_xyz, S_SAT[j]->P);
		}

		for (int j = 0; j < sat_num; j++)
		{
			double sinE = sin(sat_elev[(int)(S_SAT[j]->PRN - 1)]);
			Cl->data[j][j] = 1;//OBS_SIGMA0 * OBS_SIGMA0 / sinE;

			L->data[j][0] = S_OBS[j]->C1 - S_STA[j]->C1 + S2[j] - S[j] - solution[3];

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

		SLS(L, A, Cl, δ, Q);

		for (int j = 0; j < 3; j++) {
			solution[j] += δ->data[j][0];
			Q->data[j][j] += SLSQ_SYS_NOI;
		}

		free_mat(A);
		free_mat(L);
		free_mat(Cl);
		L = malloc_mat(sat_num, 1);
		A = malloc_mat(sat_num, 4);
		Cl = malloc_mat(sat_num, sat_num);

	}
	free_mat(L); free_mat(A); free_mat(Cl);
	free_mat(r); free_mat(δ);
	free_mat(Cx);
	return true;
}
FILE * ccsl;
FILE * ccsl2;
FILE * outfp;

void output_1()
{
	for (int i = 0; i < GPS_SAT_NUM; i++)
	{
		fprintf(ccsl, "%lf\t", phase_rate_value_1[i]);
		fprintf(ccsl2, "%lf\t", phase_rate_value_2[i]);
	}
	fprintf(ccsl, "\n");
	fprintf(ccsl2, "\n");
}

void output_2()
{
	fprintf(outfp, "%lf\t%lf\t%lf\n", solution[0], solution[1], solution[2]);
}

void check_observation()
{
	ccsl = fopen("ccsl.txt", "w");
	ccsl2 = fopen("ccsl2.txt", "w");
	ofp = fopen(OBS_FILE, "rb");
	nfp = fopen(SAT_FILE, "rb");
	sfp = fopen(STA_FILE, "rb");

	while (!feof(ofp) && !feof(nfp) && !feof(sfp))
	{
		if (!fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp)) break;
		if (!fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp)) break;
		if (!fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp)) break;
		inverse_phase(OBS);
		inverse_phase(STA);
		overall_check();
		append_history(ev);

		output_1();

		overall_reset();
	}

	_fcloseall();
}

void solve()
{
	ofp = fopen(OBS_FILE, "rb");
	nfp = fopen(SAT_FILE, "rb");
	sfp = fopen(STA_FILE, "rb");
	if (METHOD == EKF)
		ekf_create(EKF_INTERVAL);

	while (!feof(ofp) && !feof(nfp) && !feof(sfp))
	{
		fread(OBS, CHANNEL_NUM, sizeof(obs_epoch), ofp);
		fread(SAT, CHANNEL_NUM, sizeof(sat_epoch), nfp);
		fread(STA, CHANNEL_NUM, sizeof(obs_epoch), sfp);

		inverse_phase(OBS);
		inverse_phase(STA);
		overall_check();

		if (OBS->Time >= sigment.start_time && OBS->Time <= sigment.last_for + sigment.start_time)
		{
			fetch_ref_sat();
			if (METHOD == SLSQ)
				slsq_solve();
			else if (METHOD == EKF)
				ekf_solve();

			output_2();
		}
		overall_reset();
	}

	_fcloseall();
}

int main()
{
	check_observation();

	outfp = fopen("out.txt", "w");
	find_perfect_sigment();

	solve();
}