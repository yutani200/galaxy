void calc_force(int n,Istar *star);
double get_dt(int n,double T,double pre_dt,double output_time,Istar *star);
void leap_frog(int n,double dt,Istar *star);
void run(int NMAX,double Tend,double output_dt,Istar *star);
