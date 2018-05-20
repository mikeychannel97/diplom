/***************************************************************************

                                XKY_HBO_4_6.h
                             -------------------
 
 ***************************************************************************/

/*
	XKY.НВО.4.6. Базовые функции

1	FAARM				Расчет фундаментальных аргументов
2	MNUT				Расчет матрицы нутации
3	PREC2000			Расчет матрицы прецессии
4	VOSX_E				Расчет времени восхода и захода КА
5	VEKDAT_E			Определение вектора параметров орбиты в апогее на заданную дату с учетом коррекции
6	Fp_VZ				Функция-параметр расчета угла места
7	SOLM				Расчет координат Солнца
8	LUNM				Расчет координат Луны
9	t_TIME				Перевод времени из сек. в системное
10	TIME_t				Перевод времени из системного в сек.
11	OB_MK				Объединение действующего и планируемого массивов коррекций 	06	оформ-ление		
12  tabk	            Чтение констант из текстового файла TKOC_2000
13  tabk_mkc	        Чтение констант из текстового файла MKC
14  Fprc	            Формирование в файле регистрации сообщения о прохождении процесса расчета комплексной программы
15  FprcEnd	            Формирование в файле регистрации завершающего сообщения об успешном расчете комплексной программы
16  FprcCancel	        Формирование в файле регистрации завершающего сообщения о невозможности расчета комплексной программы
17  get_ibm_sys_times70	Вычисление текущей даты  и системного времени

*/

extern void FAARM (double MJD, double *L_b, double *L, double *L_1, double *F_b, double *D_b);
extern void MNUT (double MJD, double Nut[3][3]);
extern void PREC2000 (double MJD, double Prec[3][3]);
extern int VEKDAT_E(KU_TimeDATA *t0, double  P0[6], KU_DateDATA *D, int tip, int Pv,
					KU_MKOR *MK, KU_TimeDATA *ta,double  Pa[6]);
extern int VEKDAT_RU(KU_TimeDATA *t0, double  P0[6], KU_DateDATA *D, double taB, int tip,
					 int Pv, KU_MKOR *MK, KU_TimeDATA *tru,double  Pru[6]);
extern void  SOLM(double mjd,  int psk, double xs[4]);
extern void LUNM(double mjd,int psk, double xl[4]);
extern KU_TimeDATA t_TIME(double t);
extern double TIME_t(KU_TimeDATA tc);
void Fprc(char str[15], char str1[15]);
void FprcEnd(char str[15], char str1[15]);
void FprcCancel(char str[15], char str1[15]);
KU_TimeDATA get_ibm_sys_times70(KU_DateDATA *D);
