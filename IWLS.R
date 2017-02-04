## Iteratively Reweighted Least Squares

fit.glm = glm(am ~ mpg, data = mtcars, family = binomial)
summary(fit.glm)

X = model.matrix(fit.glm)
y = mtcars$am

r = 100
beta.mat = matrix(NA, ncol = 2, nrow = r)

beta.mat[1, ] = c(0, 0)

for (i in 2:r) {
  mu.vec = 1 / (1 + exp(-(X %*% beta.mat[i - 1, ])))
  sigma.mat = diag(c(mu.vec * (1 - mu.vec)))
  nu.vec = X %*% beta.mat[i - 1, ] + (y - mu.vec) / diag(sigma.mat)
  beta.mat[i, ] = solve(t(X) %*% (sigma.mat) %*% X) %*% t(X) %*% sigma.mat %*% nu.vec
}

beta.mat[r, ]
sqrt(diag(solve(t(X) %*% sigma.mat %*% X)))
summary(fit.glm)

par(mfrow = c(2, 1))
plot(beta.mat[1:10, 1], type = "l")
plot(beta.mat[1:10, 2], type = "l")
par(mfrow = c(1, 1))
  
