#include <assert.h>
#include <pthread.h>

int x = 0;

void *t1(void* arg)
{
}

void *t2(void* arg)
{
	x = 2;
}

void *t3(void* arg)
{
	x = 3;
}

int
main(void)
{
  pthread_t id1, id2, id3;

  pthread_create(&id1, NULL, t1, NULL);
  pthread_create(&id2, NULL, t2, NULL);
  pthread_create(&id3, NULL, t3, NULL);

  pthread_join(id1, NULL);
  pthread_join(id2, NULL);

  assert(x != 3);

  return 0;
}
