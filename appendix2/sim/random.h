double random_double (void);
unsigned long seed_rng (void);  /* pick "random" seed based on /dev/urandom or time */
void set_rng_seed(unsigned long newseed);
