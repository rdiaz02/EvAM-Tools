
/*	
    queue.h - integer queue, FIFO
*/


#define TRUE    1
#define FALSE   0

//typedef int bool;



#define MAX_QUEUESIZE 10000

typedef struct {
    int q[MAX_QUEUESIZE+1];
    int head;
    int tail;
    int size;
} queue;


void init_queue(queue *q);
void enqueue(queue *q, int x);
int  dequeue(queue *q);
int  empty(queue *q);








