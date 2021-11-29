//compile with: g++ -lpthread <sourcename> -o <executablename>

// Implementation of priority ceiling protocol

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/types.h>

//code of periodic tasks
void task1_code( );
void task2_code( );
void task3_code( );
void task4_code( );



//characteristic function of the thread, only for timing and synchronization
//periodic tasks
void *task1( void *);
void *task2( void *);
void *task3( void *);
void *task4( void *);



// initialization of mutexes and conditions (only for aperiodic scheduling)

#define INNERLOOP 100
#define OUTERLOOP 2000
#define NPERIODICTASKS 4
#define NTASKS NPERIODICTASKS 


struct sched_param parameters[NTASKS];
struct timespec next_arrival_time[NTASKS];
long int periods[NTASKS];
double WCET[NTASKS];
int missed_deadlines[NTASKS];
pthread_attr_t attributes[NTASKS];
pthread_t thread_id[NTASKS];
// global variables to be written and read inside critical sections
int t1t2 = 0;
int t2t3 = 0;
int t1t4 = 0;

// variables for priority ceiling
pthread_mutex_t mutext1t2 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutext1t4 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutext2t3 = PTHREAD_MUTEX_INITIALIZER;

//Bij:j-th Critical section of i-th task  
double B11;
double B12;
double B21;
double B22;
double B31;
double B41;

// Bi :the longest Critical sections that block Ti, directly or indirectly
double B1;
double B2;
double B3;
double B4;
//B4[] this is empty because T4 is the one with the lowest priority


int main()
{	

  	//you can already order them according to their priority; 
	//if not, you will need to sort them
  	periods[0]= 80000000; //in nanoseconds
  	periods[1]= 100000000; //in nanoseconds
  	periods[2]= 160000000; //in nanoseconds
  	periods[3]= 200000000; //in nanoseconds

  	
	//this is not strictly necessary, but it is convenient to
	//assign a name to the maximum and the minimum priotity in the
	//system. We call them priomin and priomax.

  	struct sched_param priomax;
  	priomax.sched_priority=sched_get_priority_max(SCHED_FIFO);
  	struct sched_param priomin;
  	priomin.sched_priority=sched_get_priority_min(SCHED_FIFO);

	// set the maximum priority to the current thread (you are required to be
  	// superuser). Check that the main thread is executed with superuser privileges
	// before doing anything else.

  	if (getuid() == 0)
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax);
	// execute all tasks in standalone modality in order to measure execution times
  	// (use gettimeofday). Use the computed values to update the worst case execution
  	// time of each task.

 	
 	// initialize time_i required to read the clock
	
  	struct timespec time_x, time_y;
	int i;
  	for (i =0; i < NPERIODICTASKS; i++)
    	{

	//we should execute each task more than one for computing the WCET
		//periodic tasks
 	     	clock_gettime(CLOCK_REALTIME, &time_x);
 	     	if (i== 0)
			task1_code();
			
		if (i == 1)
			task2_code();

      		if (i == 2)
			task3_code();
			
		if (i == 3)
			task4_code();
      		
      		clock_gettime(CLOCK_REALTIME, &time_y);
      		
      		// Below (inside the execution of task-i functions) I have ricomputed WCET for a more precision value. Every 10 iterations w.r.t 100 in total
		
		
		WCET[i]= 1000000000*(time_y.tv_sec - time_x.tv_sec)
			       +(time_y.tv_nsec-time_x.tv_nsec);
		
      		printf("\nWorst Case Execution Time %d=%f \n", 0, WCET[i]);
		
		}
	// compute the set of critical blocking sections for each task 
    		
		
		double B1Star[3] = {B21,B22,B41};
		double B2Star[2] = {B31,B41};
		double B3Star[1] = {B41};
		
    	
    	// compute the longest critical section for each set
    	double max1 = B1Star[0];
    	for(int i= 0; i< sizeof(B1Star); i++){
    		
    		if (B1Star[i] > max1){
    			max1 = B1Star[i];
    		}
    		B1 = max1;
    		
    	}
    	
    	double max2 = B2Star[0];
    	for(int i= 0; i< sizeof(B2Star); i++){
    		
    		if (B2Star[i] > max2){
    			max2 = B2Star[i];
    		}
    		B2 = max2;
    		
    	}
    	
    	double max3 = B3Star[0];
    	for(int i= 0; i< sizeof(B3Star); i++){
    		
    		if (B3Star[i] > max3){
    			max3 = B3Star[i];
    		}
    		B3 = max3;
    		
    	}
    	
    	
    	
    	// compute U 
    	
	double U1 = WCET[0]/periods[0]+ B1/periods[0];
	
	double U2 = WCET[0]/periods[0]+WCET[1]/periods[1]+ B2/periods[1];
	
	double U3 = WCET[0]/periods[0]+WCET[1]/periods[1]+WCET[2]/periods[2] + B3/periods[2];
	
	double U4 = WCET[0]/periods[0]+WCET[1]/periods[1]+WCET[2]/periods[2] + WCET[3]/periods[3] ;

    	
    	// compute Ulub by considering the fact that we have harmonic relationships between periods
	double UlubT1 = 1;
	double UlubT2 = 2*(pow(2.0,(1.0/2)) -1);
	double UlubT3 = 3*(pow(2.0,(1.0/3)) -1);
	double UlubT4 = 4*(pow(2.0,(1.0/4)) -1);
    	
	//if there are no harmonic relationships, use the following formula instead
	//double Ulub = NPERIODICTASKS*(pow(2.0,(1.0/NPERIODICTASKS)) -1);
	
	//check the sufficient conditions: if they are not satisfied, exit  
  	if (U1 > UlubT1)
    	{
      		printf("\n U1=%lf UlubT1=%lf Non schedulable Task Set", U1, UlubT1);
      		return(-1);
    	}
  	printf("\n U1=%lf UlubT1=%lf Scheduable Task Set", U1, UlubT1);
  	fflush(stdout);
  	sleep(5);
  	
  	
  	if (U2 > UlubT2)
    	{
      		printf("\n U2=%lf UlubT2=%lf Non schedulable Task Set", U2, UlubT2);
      		return(-1);
    	}
  	printf("\n U2=%lf UlubT2=%lf Scheduable Task Set", U2, UlubT2);
  	fflush(stdout);
  	sleep(5);
  	
  	
  	if (U3 > UlubT3)
    	{
      		printf("\n U3=%lf UlubT3=%lf Non schedulable Task Set", U3, UlubT3);
      		return(-1);
    	}
  	printf("\n U3=%lf UlubT3=%lf Scheduable Task Set", U3, UlubT3);
  	fflush(stdout);
  	sleep(5);
  	
  	
  	if (U4 > UlubT4)
    	{
      		printf("\n U4=%lf UlubT4=%lf Non schedulable Task Set \n", U4, UlubT4);
      		return(-1);
    	}
  	printf("\n U4=%lf UlubT4=%lf Scheduable Task Set \n", U4, UlubT4);
  	fflush(stdout);
  	sleep(5);

	// set the minimum priority to the current thread: this is now required because 
	//we will assign higher priorities to periodic threads to be soon created
	//pthread_setschedparam
	
	
  	// set the attributes of each task, including scheduling policy and priority
  	
	if (getuid() == 0)
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);
	int j;
	for (j =0; j < NPERIODICTASKS; j++)
    	{
    	
    	//initialize the attribute structure of task i
      	pthread_attr_init(&(attributes[j]));

//set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
//not inherited from the main thread (pthread_attr_setinheritsched) 
      	pthread_attr_setinheritsched(&(attributes[j]), PTHREAD_EXPLICIT_SCHED);
      
// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
	pthread_attr_setschedpolicy(&(attributes[j]), SCHED_FIFO);

//properly set the parameters to assign the priority inversely proportional to the period
      	//parameters[j].sched_priority = priomin.sched_priority+NTASKS - j;
parameters[j].sched_priority = sched_get_priority_max(SCHED_FIFO) -j;
//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)

      	pthread_attr_setschedparam(&(attributes[j]), &(parameters[j]));
      	}
      	
      	pthread_mutexattr_t mymutexattr1; 
      	pthread_mutexattr_t mymutexattr2; 
      	pthread_mutexattr_t mymutexattr3; 
      	 
      	pthread_mutexattr_init(&mymutexattr1);	
      	pthread_mutexattr_init(&mymutexattr2);	
      	pthread_mutexattr_init(&mymutexattr3);	
      	
    	pthread_mutexattr_setprotocol(&mymutexattr1, PTHREAD_PRIO_PROTECT);
    	pthread_mutexattr_setprotocol(&mymutexattr2, PTHREAD_PRIO_PROTECT);
    	pthread_mutexattr_setprotocol(&mymutexattr3, PTHREAD_PRIO_PROTECT);
    	 
    	pthread_mutexattr_setprioceiling(&mymutexattr1,parameters[0].sched_priority); 
    	
    	pthread_mutexattr_setprioceiling(&mymutexattr2,parameters[0].sched_priority); 
    	
    	pthread_mutexattr_setprioceiling(&mymutexattr3,parameters[1].sched_priority); 
      	
    		
      	
	pthread_mutex_init(&mutext1t2, &mymutexattr1); 
	
	pthread_mutex_init(&mutext1t4, &mymutexattr2); 
	
	pthread_mutex_init(&mutext2t3, &mymutexattr3); 
  	
      		
	
  	//declare the variable to contain the return values of pthread_create	
  	int iret[NTASKS];

	//declare variables to read the current time
	struct timespec time_f;
	clock_gettime(CLOCK_REALTIME, &time_f);

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NPERIODICTASKS; i++)
    	{
		long int next_arrival_nanoseconds = time_f.tv_nsec + periods[i];
		//then we compute the end of the first period and beginning of the next one
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_f.tv_sec + next_arrival_nanoseconds/1000000000;
       		missed_deadlines[i] = 0;
    	}

	

	// create all threads(pthread_create)
  	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
  	iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);

  	// join all threads (pthread_join)
  	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
  	pthread_join( thread_id[3], NULL);

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NTASKS; i++)
    	{
      		printf ("\nMissed Deadlines Task %d= %d", i, missed_deadlines[i], "\n");
		fflush(stdout);
    	}
    	pthread_mutexattr_destroy(&mymutexattr1);
    	pthread_mutexattr_destroy(&mymutexattr2);
    	pthread_mutexattr_destroy(&mymutexattr3);
  	exit(0);
  	
}

// application specific task_1 code
void task1_code()
{
	int i,j;
	double uno;
	struct timespec time_1, time_2;
	struct timespec time_3, time_4;
	
	
	
	//this double loop with random computation is only required to waste time
	
	//print the id of the current task
  	printf(" 1[ "); fflush(stdout);
	
	pthread_mutex_lock(&mutext1t2);
	clock_gettime(CLOCK_REALTIME, &time_1);
  	for (i = 0; i < OUTERLOOP; i++)
    	{
    		for (j = 0; j < INNERLOOP; j++)
		{
			uno = rand()*rand()%10;
		}
	}
			
	t1t2++;
	clock_gettime(CLOCK_REALTIME, &time_2);
	pthread_mutex_unlock(&mutext1t2);
			
	// I computeed it even if it is not usefull	
	B11= 1000000000*(time_2.tv_sec - time_1.tv_sec)
			       +(time_2.tv_nsec-time_1.tv_nsec);
			
	
		
    	pthread_mutex_lock(&mutext1t4);
	clock_gettime(CLOCK_REALTIME, &time_3);
  	for (i = 0; i < OUTERLOOP; i++)
    	{
    		//print the id of the current task
  	

      		for (j = 0; j < INNERLOOP; j++)
		{
		uno = rand()*rand()%10;
		}
	
  	}
  	t1t4++;
  	//print the id of the current task
  	printf(" ]1 "); fflush(stdout);
			
	clock_gettime(CLOCK_REALTIME, &time_4);
    	pthread_mutex_unlock(&mutext1t4);
    	
    	// I computeed it even if it is not usefull	
    	B12= 1000000000*(time_4.tv_sec - time_3.tv_sec)
			       +(time_4.tv_nsec-time_3.tv_nsec);	
    		
	
  	
}

//thread code for task_1 (used only for temporization)
void *task1( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

   	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	struct timespec time_a, time_b;
	double temp;
  	for (i=0; i < 100; i++)
    	{
      		// execute application specific code
      		// every 10 iterations I recalculate the WCET for a better precision in the computation
      	
      	if (i%10 == 0){
      		clock_gettime(CLOCK_REALTIME, &time_a);	
			task1_code();
		clock_gettime(CLOCK_REALTIME, &time_b);
		temp = 1000000000*(time_b.tv_sec - time_a.tv_sec)+(time_b.tv_nsec-time_a.tv_nsec);
			
		if (temp > WCET[0]){
			WCET[0] = temp;
			}
	}
	else
	
		task1_code();
		
		// sleep until the end of the current period (which is also the start of the
		// new one
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL);

		// the thread is ready and can compute the end of the current period for
		// the next iteration
 		
		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task2_code()
{
	struct timespec time_5, time_6;
	struct timespec time_7, time_8;
	int i,j;
	int choice1;
	double uno;
	
	
	//print the id of the current task
  	printf(" 2[ "); fflush(stdout);
	pthread_mutex_lock(&mutext1t2);
	clock_gettime(CLOCK_REALTIME, &time_5);
  	for (i = 0; i < OUTERLOOP; i++)
    	{
    	
      		for (j = 0; j < INNERLOOP; j++)
		{
		uno = rand()*rand()%10;
		}
		
		
	}
	if (t1t2 >= 0){
		choice1 = 1;
	}else{
		choice1 = 0;
	}
	clock_gettime(CLOCK_REALTIME, &time_6);
			pthread_mutex_unlock(&mutext1t2);
	B21= 1000000000*(time_6.tv_sec - time_5.tv_sec)
			       +(time_6.tv_nsec-time_5.tv_nsec);
	
	
	pthread_mutex_lock(&mutext2t3);
	clock_gettime(CLOCK_REALTIME, &time_7);
    	for (i = 0; i < OUTERLOOP; i++)
    	{
    	
      		for (j = 0; j < INNERLOOP; j++)
		{
		uno = rand()*rand()%10;
		}
	
    	}
	
    	t2t3++;		
	clock_gettime(CLOCK_REALTIME, &time_8);
	pthread_mutex_unlock(&mutext2t3);
	
	//print the id of the current task
  	printf(" ]2 "); fflush(stdout);		
			
	B22= 1000000000*(time_8.tv_sec - time_7.tv_sec)
			       +(time_8.tv_nsec-time_7.tv_nsec);
		
	
}
	

void *task2( void *ptr )
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
	struct timespec time_c, time_d;
	double temp;
	
  	for (i=0; i < 100; i++)
    	{
      		
      		if (i%10 == 0){
      		clock_gettime(CLOCK_REALTIME, &time_c);	
		task2_code();
		clock_gettime(CLOCK_REALTIME, &time_d);
		temp = 1000000000*(time_d.tv_sec - time_c.tv_sec)+(time_d.tv_nsec-time_c.tv_nsec);
			
		if (temp > WCET[1]){
			WCET[1] = temp;
			}
	}
	else
		task2_code();
		
	clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task3_code()
{
	struct timespec time_9, time_10;
	int i,j;
	int choice2;
	double uno;
	
	
	
	//print the id of the current task
  	printf(" 3[ "); fflush(stdout);
	pthread_mutex_lock(&mutext2t3);
      	clock_gettime(CLOCK_REALTIME, &time_9);		
  	for (i = 0; i < OUTERLOOP; i++)
    	{
    		for (j = 0; j < INNERLOOP; j++){
      		uno = rand()*rand()%10;	
		}
	}
	
	if (t2t3 >= 0 ){
		choice2 = 1;
	} 
	else{
		choice2 = 0;
	}
	clock_gettime(CLOCK_REALTIME, &time_10);
    	pthread_mutex_unlock(&mutext2t3);
    	
    	//print the id of the current task
  	printf(" ]3 "); fflush(stdout);		
    	B31= 1000000000*(time_10.tv_sec - time_9.tv_sec)
			       +(time_10.tv_nsec-time_9.tv_nsec);
    			
}

void *task3( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
	struct timespec time_e, time_f;
	double temp;
  	for (i=0; i < 100; i++)
    	{
      	
      	if (i%10 == 0){
      		clock_gettime(CLOCK_REALTIME, &time_e);	
		task3_code();
		clock_gettime(CLOCK_REALTIME, &time_f);
		temp = 1000000000*(time_f.tv_sec - time_e.tv_sec)+(time_f.tv_nsec-time_e.tv_nsec);
			
		if (temp > WCET[2]){
			WCET[2] = temp;
			}
	}
	else
	
		task3_code();
		

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
		next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}



void task4_code()
{
	struct timespec time_11, time_12;
	int i,j;
	int choice3;
	double uno;
	
	
	//print the id of the current task
  	printf(" 4[ "); fflush(stdout);
	
	pthread_mutex_lock(&mutext1t4);
      	clock_gettime(CLOCK_REALTIME, &time_11);			
  	for (i = 0; i < OUTERLOOP; i++)
    	{
    		
      		for (j = 0; j < INNERLOOP; j++){
      		uno = rand()*rand()%10;
      			
		}
	
	}
	if(t1t4>= 0){
		choice3 = 1;
	}else{
		choice3 = 0;
	}
	clock_gettime(CLOCK_REALTIME, &time_12);	
    	pthread_mutex_unlock(&mutext1t4);
    	
    	//print the id of the current task
  	printf(" ]4 "); fflush(stdout);		
    	
    	B41= 1000000000*(time_12.tv_sec - time_11.tv_sec)
			       +(time_12.tv_nsec-time_11.tv_nsec);
    	
	
}

void *task4( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
	struct timespec time_g, time_h;
	double temp;
  	for (i=0; i < 150; i++)
    	{
    	
      		if (i%10 == 0){
      		clock_gettime(CLOCK_REALTIME, &time_g);	
		task4_code();
		clock_gettime(CLOCK_REALTIME, &time_h);
		temp = 1000000000*(time_h.tv_sec - time_g.tv_sec)+(time_h.tv_nsec-time_g.tv_nsec);
			
		if (temp > WCET[3]){
			WCET[3] = temp;
			}
	}
	else
	
		task4_code();
		

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
		next_arrival_time[3].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec= next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}



