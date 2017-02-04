## Coronary Heart Disease (CHD) - Hosmer and Lemeshow (1989)
##  100 patients, by age and CHD status

chd.df=data.frame(age=c(20,23,24,25,25,26,26,28,28,29,30,30,30,30,30,30,32,32,33,33,34,34,34,34,34,35,35,36,36,36,37,37,37,38,38,39,39,40,40,41,41,42,42,42,42,43,43,43,44,44,44,44,45,45,46,46,47,47,47,48,48,48,49,49,49,50,50,51,52,52,53,53,54,55,55,55,56,56,56,57,57,57,57,57,57,58,58,58,59,59,60,60,61,62,62,63,64,64,65,69),
                  chd=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,1,0,1,0,1,0,0,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,1))


## 95% prior probability proportion of chd is between 0.1 and 0.5
b0 = 1 / 2 * (log(0.1 / (1 - 0.1)) + log(0.5 / (1 - 0.5)))
s0 = 1 / (2 * 1.960) * (log(0.5 / (1 - 0.5)) - log(0.1 / (1 - 0.1)))

## Age has s = 11.4
## Increasing age by 10 years equally likely to multiply odds of chd by a
##  a factor above 2 or below 2. With 95% prior probability it will be by no
##  more than 8 though.
b.age = log(2) / 10
s.age = 1 / (1.645 * 10) * log(8 / 2)

b0.vec = c(b0, b.age)
V0.mat = diag(c(s0, s.age)^2)

fit.ml = glm(chd ~ scale(age, scale = FALSE), data = chd.df, family = binomial)
summary(fit.ml)

beta.ml = coef(fit.ml)
V.ml = vcov(fit.ml)

V1.mat = solve(solve(V0.mat) + solve(V.ml))
b1.vec = V1.mat %*% solve(V0.mat) %*% b0.vec + V1.mat %*% solve(V.ml) %*% beta.ml

L = chol(V1.mat)

nsteps = 5000
beta.mat = matrix(NA, ncol = 2, nrow = nsteps)
beta.mat[1, ] = b1.vec + L %*% rt(2, 4)

y = chd.df$chd
X = model.matrix(fit.ml)

library(mvtnorm)

for (i in 2:nsteps) {
  betaprime.vec = c(b1.vec + L %*% rt(2, 4))
  p.prime = 1 / (1 + exp(-(X %*% betaprime.vec)))
  loglik.prime = sum(dbinom(y, 1, p.prime, log = TRUE))
  prior.prime = sum(dnorm(betaprime.vec, b0.vec, sqrt(diag(V0.mat)), log = TRUE))
  candidate.prime = dmvt(betaprime.vec, b1.vec, V1.mat, 4, log = TRUE)
#  candidate.prime = sum(dt((betaprime.vec - b1.vec) / sqrt(diag(V1.mat)), 4, log = TRUE))

  p.current = 1 / (1 + exp(-(X %*% beta.mat[i - 1, ])))
  loglik.current = sum(dbinom(y, 1, p.current, log = TRUE))
  prior.current = sum(dnorm(beta.mat[i - 1, ], b0.vec, sqrt(diag(V0.mat)), log = TRUE))
  candidate.current = dmvt(beta.mat[i - 1, ], b1.vec, V1.mat, 4, log = TRUE)
#  candidate.current = sum(dt((beta.mat[i - 1, ] - b1.vec) / sqrt(diag(V1.mat)), 4, log = TRUE))

  alpha = exp(loglik.prime + prior.prime + candidate.current - (loglik.current + prior.current + candidate.prime))

  if (runif(1) <= alpha) beta.mat[i, ] = betaprime.vec
  else beta.mat[i, ] = beta.mat[i - 1, ]
}

plot(beta.mat[1:200, 1], type = "l")
acf(beta.mat[-(1:200), 1])

plot(beta.mat[1:200, 2], type = "l")
acf(beta.mat[-(1:200), 2])

mean(beta.mat[-(1:10), 1][(1:(nsteps - 10)) %% 5 == 0])
sd(beta.mat[-(1:10), 1][(1:(nsteps - 10)) %% 5 == 0])
hist(beta.mat[-(1:10), 1][(1:(nsteps-10)) %% 5 == 0], freq = FALSE)
lines(seq(-2, 1, 0.01), dnorm(seq(-2, 1, 0.01), b1.vec[1], sqrt(V1.mat[1, 1])))

mean(beta.mat[-(1:10), 2][(1:(nsteps-10)) %% 5 == 0])
sd(beta.mat[-(1:10), 2][(1:(nsteps - 10)) %% 5 == 0])
hist(beta.mat[-(1:10), 2][(1:(nsteps - 10))%%5 == 0], freq = FALSE)
lines(seq(0, 0.3, 0.001),dnorm(seq(0, 0.3, 0.001), b1.vec[2], sqrt(V1.mat[2, 2])))





