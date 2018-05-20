#include "stdafx.h"
/***************************************************************************

                                XKY_HBO_4_2.h
                             -------------------
 
 ***************************************************************************/

/*

	XKY.НВО.4.2.Прогнозирование орбитального движения КА

1	INT2000		Прогнозирование орбитального движения КА на основе
					численного интегрирования уравнений движения	
2	INTUM		Прогноз орбитального движения до заданного момента времени
3	INTUUM		Прогноз орбитального движения до заданного аргумента широты
4	INTUS		Прогнозирование движения КА с учетом управляющих ускорений
5	INTUKOR		Прогноз орбитального движения с учетом управляющих
					ускорений до заданного времени
6	INTUUKOR	Прогноз орбитального движения с учетом управляющих
					ускорений до заданного аргумента широты
7	INTK		Определение движения КА с учетом протяженных коррекций
8	INTK_1		Определение движения КА с учетом одной протяженной коррекции
9	INTKU		Прогнозирование движения КА до заданного аргумента широты
					с учетом управляющих воздействий
10	APSIDM		Прогноз до апсидальной точки
11	APSIDK		Прогноз движения КА в апсидальную точку с учетом
					управляющих воздействий
12	INJ2000		Прогноз орбитального движения КА в системе координат J2000				
13	INJUS		Прогноз орбитального движения КА в системе координат J2000
					с учётом управляющих воздействий				
14	INJUM		Прогноз орбитального движения КА до заданного момента
					времени в системе координат J2000
15	INJUUM		Прогноз орбитального движения КА до заданного аргумента
					широты в системе координат J2000
16	INJUKOR		Прогноз орбитального движения КА до заданного момента
					времени в системе координат J2000 c учетом управляющих ускорений
17	INJUUKOR	Прогноз орбитального движения КА до заданного аргумента
					широты в системе координат J2000 c учетом управляющих ускорений
18	DRAKM		Расчет драконического периода
19	PROV		Прогноз орбитального движения до заданного времени с учетом
					действующих и планируемых коррекций

*/

int INT2000		(KU_TimeDATA *t_n,TimeDATA_u *t_k, double p[], int dot,int priz,int tip);
int INTUM		(KU_TimeDATA *tn, KU_TimeDATA *tk, double pn[], double pk[],int s);
int INTUUM		(KU_TimeDATA *tn, KU_TimeDATA *tk, int kvit, double uk, double pn[],double pk[],int s);
int INTUS		(KU_TimeDATA *t_n,TimeDATA_u *t_k, double p[], int dot,int priz,int tip,double stw[3]);
int INTUKOR		(KU_TimeDATA *tn, KU_TimeDATA *tk, double pn[], double pk[],int s,double stw[3]);
int INTUUKOR	(KU_TimeDATA *tn, KU_TimeDATA *tk, int kvit, double uk,double pn[], double pk[],int s,double stw[3]);
void INTK(KU_TimeDATA *tn, KU_TimeDATA *tk, double Pn[6], int tip, KU_MKOR *ST_KOR, int *kik, double Pk[6]);
void INTK_1(KU_TimeDATA *tn, KU_TimeDATA *tk, double Pn[6],int tip,
			KU_TimeDATA *tndk, KU_TimeDATA *tkdk, double W[3], double Pk[6]);
int INTKU(KU_TimeDATA *t ,double p[6],int kv,double uk,int tip, KU_TimeDATA *tk,double pk[6],KU_MKOR *STK);
int APSIDM(	KU_TimeDATA* tn, double pn[], int tip, int j, int ap, KU_TimeDATA *tk, double pk[] );
int APSIDK(	KU_TimeDATA* tn, double pn[], int tip, int j, int ap, KU_MKOR *STK, KU_TimeDATA *tk, double pk[] );
int DRAKM(KU_TimeDATA t,double p[6],int tip,double *tdr);
void PROV(KU_TimeDATA todn,double podn[6],int tip,KU_TimeDATA t,double p[6],
	        KU_MKOR  MDKOR,KU_MKOR  MPKOR);
void RKSHM(KU_TimeDATA *t, int *nvit, double pp[], double uscu[],
	       int sbr, int dot, int tip, int priz,  double pred[]);
void SECUSM(KU_TimeDATA *t, double tk, int *data, double pp[],  double uscu[], int tip,
	        int priz, double pred[]);
void TUZM(double u,double *sec,double pred[]);
void RKSHUS(KU_TimeDATA *t, int *nvit, double pp[], double uscu[],
	           int sbr, int dot, int tip, int priz,  double pred[]);
void SECUSS(KU_TimeDATA *t, double tk, int *data, double pp[],  double uscu[], int tip,
	            int priz, double pred[]);
void PRAVM(int data,double rr[],double z[],double uscu[], int tip, int   priz,double pred[]);
void PRAVUS(int data,double rr[],double z[],double uscu[], int tip, int   priz,double pred[]);
void SDSTW(KU_TimeDATA *ts, double fr, double *s, double *t, double *w,
		    double stw[]);
void LS2000(KU_TimeDATA *ts, int l, double *s, double *t, double *w, double stw[]);
void GEOM(double sek, double *gs, double *gt, double *gw, double sm0,double stw[]);




