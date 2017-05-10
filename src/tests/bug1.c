#include <assert.h>
#include <pthread.h>

int x = 0;
pthread_mutex_t m;

void *t2(void* arg)
{
	pthread_mutex_lock(&m);
	x = 2;
	pthread_mutex_unlock(&m);
}

void *t3(void* arg)
{
	pthread_mutex_lock(&m);
//	x = 3;
	pthread_mutex_unlock(&m);
}

int
main(void)
{
  pthread_t id1, id2, id3;
  pthread_mutex_init(&m, NULL);

  pthread_create(&id2, NULL, t2, NULL);
  pthread_create(&id3, NULL, t3, NULL);

  pthread_mutex_lock(&m);
  x = 4;
  assert(x == 4);
  pthread_mutex_unlock(&m);

  return 0;
}
