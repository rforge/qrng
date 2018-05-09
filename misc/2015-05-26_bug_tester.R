require(qrng)

sobol(10, d=5, randomize=TRUE)
u <- sobol(1e7, d=4)

for(i in 1:10)
    u <- sobol(1e7, d=4)
