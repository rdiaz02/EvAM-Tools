
/*	
    queue.c - integer queue, FIFO
*/

#include "queue.h"

#include <stdio.h>
#include <stdlib.h>


void init_queue(queue *q)
{
    q->head = 0;
    q->tail = MAX_QUEUESIZE-1;
    q->size = 0;
}


void enqueue(queue *q, int x)
{
    if (q->size >= MAX_QUEUESIZE)
    {
	fprintf(stderr, "Error: queue overflow at enqueue  x = %d\n", x);
	exit(1);
    }
    else 
    {
	q->tail = (q->tail+1) % MAX_QUEUESIZE;
	q->q[q->tail] = x;    
	q->size = q->size + 1;
    }
}


int dequeue(queue *q)
{
  int x = 0;
    
    if (q->size <= 0) 
	fprintf(stderr, "Warning: empty queue dequeue.\n");
    else 
    {
	x = q->q[q->head];
	q->head = (q->head+1) % MAX_QUEUESIZE;
	(q->size)--;
    }
    
    return(x);
}


int empty(queue *q)
{
    return (q->size <= 0);
}

