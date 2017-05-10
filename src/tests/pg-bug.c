/* Adapted from PGSQL benchmark from http://link.springer.com/chapter/10.1007%2F978-3-642-37036-6_28 */

#include <stdbool.h>
#include <assert.h>
#include <pthread.h>


volatile _Bool latch1 = true;
volatile _Bool latch2 = false;

void* worker_1(void* arg)
{
    latch2 = true;

    if (!latch1);
    if (!latch1);
    latch2 = true;
  return NULL;
}

void* worker_2(void* arg)
{
    if (!latch2);
    latch2 = false;
    latch1 = true;

    if (!latch2);
    latch2 = false;
    latch1 = true;
  return NULL;
}

int main() {
  pthread_t t1, t2;

  pthread_create(&t1, 0, worker_1, NULL);
  pthread_create(&t2, 0, worker_2, NULL);

  return 0;
}
