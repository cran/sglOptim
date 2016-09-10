library(sglOptim)

# warnings = errors
options(warn=2)

data(TestData)
x <- test.data$x
y <- test.data$y
grp <- test.data$grp

weights <- rep(1/nrow(x), nrow(x))
sampleGrouping <- grp
covariateGrouping <- factor(1:ncol(x))
groupWeights <- c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x))
alpha <- 0
d <- 50L
lambda.min <- 2
algorithm.config <- sgl.standard.config

# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
lambda <- sgl_lambda_sequence("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
fit1a <- sgl_fit("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)

# Predict test
res1a <- sgl_predict("sgl_test_dense", "sglOptim", fit1a, data)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
fit1b <- sgl_fit("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)

# Predict test
res1b <- sgl_predict("sgl_test_sparse", "sglOptim", fit1b, data)


if(max(abs(fit1a$beta[[25]]-fit1b$beta[[25]])) > 1e-5) stop()
if(max(abs(res1a$responses$link[[25]]-res1b$responses$link[[25]])) > 1e-5) stop()

## Navigate tests

# nmod
if(nmod(fit1a) != d) stop()

# models
m <- models(fit1a)
if(length(m) != d) stop()
if(any(colnames(m[[1]]) != colnames(x))) stop()
if(any(rownames(m[[1]]) != levels(grp))) stop()

# coef
if(length(coef(fit1a)) != d) stop()

# features and parameters
features(fit1a)
parameters(fit1a)
features_stat(fit1a)
parameters_stat(fit1a)

# print
sgl_print(fit1a)
sgl_print(res1a)

# test names on output
if(any(colnames(fit1a$beta[[10]]) != colnames(x))) stop()

if(any(rownames(fit1a$beta[[10]]) != levels(grp))) stop()

if(any(colnames(res1a$responses$link[[10]]) != rownames(x))) stop()
