extern void __VERIFIER_error() __attribute__ ((__noreturn__));

#include <pthread.h>

void __VERIFIER_assert(int expression) { if (!expression) { ERROR: __VERIFIER_error();}; return; }

char v;

void *thread2(void *arg)
{
  v = 'X';
}

void *thread3(void *arg)
{
  v = 'Y';
}

int main(void)
{
  pthread_t t1, t2, t3, t4, t5;

  pthread_create(&t3, 0, thread3, 0);
  pthread_create(&t5, 0, thread2, 0);

  pthread_join(t3, 0);
  pthread_join(t5, 0);

  __VERIFIER_assert(v == 'X'); // <-- wrong, the only thread that writes 'Y' can be the last to write

  return 0;
}

