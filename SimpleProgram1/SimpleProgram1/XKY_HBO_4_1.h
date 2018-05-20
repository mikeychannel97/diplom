#include "stdafx.h"
/***************************************************************************

                                XKY_HBO_4_1.h
                             -------------------
 
 ***************************************************************************/

/*

	XKY.НВО.4.1. Пересчет систем координат и времени

1	RIPM		Пересчет параметров орбиты  в специальную оскулирующую
2	PIRM		Пересчет параметров орбиты из специальной оскулирующей системы координат
3	PERUGM		Пересчет углов из экваториальной системы координат в меридиональную и обратно
4	PEING		Пересчет параметров из инерциальной системы координат в систему ПЗ-90
5	STISF		Пересчет параметров из прямоугольной стационарной системы координат в сферическую систему	
6	POSTМ		Пересчет вектора орбитальных параметров в прямоугольную стационарную  систему координат
7	OSКPAR		Расчет оскулирующих параметров КА
8	GEOCМ		Определение геоцентрических координат КА
9	GEOD		Определение геодезических координат КА
10	DIRM		Пересчет из декартовой системы координат в параметры орбиты
11	DECORM		Пересчет из параметров орбиты в декартову систему координат
12	PETJ		Пересчет вектора из текущей системы координат в систему координат J2000
13	PEJT		Пересчет вектора из системы координат J2000 в текущую систему координат
14	SLOGM		Сложение системных времен
15	SCOVM		Сравнение системных времен
16	ZVEWS		Расчет среднего звездного времени
17	ZVEWI		Расчет истинного звездного времени
18	KDS2000		Пересчет из календарного времени в системное
19	SKD2000		Пересчет системного времени 2000г в календарное
20	MJS2000		Пересчет модифицированной юлианской даты в системное время от 2000 года
21	SMJ2000		Расчет модифицированной юлианской даты по системному времени от 2000 года
22	KDMJ		Расчет модифицированной юлианской даты
23	INFOE		Расчет набора информационных параметров
24	PEPAMK		Пересчет орбитальных параметров в станционную систему координат
25	PEGIN		Пересчет параметров из системы координат ПЗ-90 в инерциальную систему
26	POSP		Пересчет оскулирующих параметров из системы координат J2000 в текущую систему координат
27	POSJ		Пересчет оскулирующих параметров из текущей системы координат систему координат J2000
28	PSORK		Пересчет из связанной системы координат в орбитальную при коррекции ОДУ
29	POSKS		Пересчет из орбитальной системы координат в связанную при коррекции ОДУ

*/


extern int RIPM(double[],double[]);
extern int PIRM(double[],double[]);
extern int PERUGM(double [],int ,double []);
extern int PEING(double Xin[6], double MJD, double DEL_UT1, double Xpz[6]);
extern int STISF(double [],double *, double *, double *, double *);
extern int POSTM(double [], double, int, double, int, double []);
extern void OSKPAR(double p[6],int tip,KU_OskPar *op);
extern void GEOCM(KU_TimeDATA t,double p[6],int tip,double DEL_UT1,
		          double*l,double*sh);
extern void GEOD(double X,double Y,double Z, double *B,double *L,double *H);
extern void DIRM(double X[6],double R[6],double Rm[6]);
extern void DECORM(double r[6],int tip, double X[6]);
extern void PETJ( KU_TimeDATA *t,double p[6],int psk,int tip,double xj[6],double pj[6]);
extern void PEJT(KU_TimeDATA *t,double pj[6],int psk,int tip,
           double x[6],double p[6]);
extern int SLOGM(KU_TimeDATA *t1,KU_TimeDATA *t2, int k, KU_TimeDATA *t3);
extern int SCOVM(KU_TimeDATA *t1,KU_TimeDATA *t2, int *a );
extern int ZVEWS(double MJD,  double DEL_UT1, double *s_m);
extern int ZVEWI(double MJD,  double DEL_UT1,  double N_alfa, double *s);
extern int KDS2000(KU_DateDATA *D, KU_TimeDATA *t);
extern int SKD2000(KU_TimeDATA *t, KU_DateDATA *D );
extern int MJS2000(double MJD, KU_TimeDATA *t);
extern int SMJ2000(KU_TimeDATA *t, double *MJD );
extern void KDMJ(KU_DateDATA *D, double *MJD);
extern int INFOE(KU_TimeDATA t,double p[6],int tip,double DEL_UT1, KU_INFOE *opi);
extern int PEPAMK(KU_TimeDATA* t,KU_TimeDATA* t_t, double rh[6],
		          int tip,  double DEL_UT1, int NPOST,  KU_MKOR *ST_KOR, double dc[3]);
extern int PEGIN(double Xpz[6], double MJD, double DEL_UT1,double Xin[6]);
extern void POSJ(KU_TimeDATA t,int v,double DEL_UT1,KU_OSKP op,
                 double *lvu,KU_OSKP *opj);
extern void POSP(KU_TimeDATA t,int v,double DEL_UT1,KU_OSKP opj,
	             double *lvu,KU_OSKP *op);
