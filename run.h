void calc_force(int n,ISTAR *star);
double get_dt(int n,double T,double pre_dt,double output_time,ISTAR *star);
void leap_frog(int n,double dt,ISTAR *star);
void initialize_run(PAP,ISTAR*);
void run(PAP, ISTAR *star);
