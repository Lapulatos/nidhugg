#include <assert.h>
#include <pthread.h>

int x = 0;

void *t1(void* arg)
{
}

void *t2(void* arg)
{
	x = 3;
	x = 4;
}

void *t3(void* arg)
{
	x = 5;
}

int
main(int argc, char **argv)
{
  pthread_t id1, id2, id3;

  pthread_create(&id1, NULL, t1, NULL);
  pthread_create(&id2, NULL, t2, NULL);
  pthread_create(&id3, NULL, t3, NULL);
  /*
  pthread_join(id1,NULL);
  pthread_join(id2,NULL);
  */

  x = 7;
  assert(x != 5);

  return 0;
}
