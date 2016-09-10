library(sglOptim)

# warnings = errors
options(warn=2)

data(TestData)
x <- test.data$x[,1:10]
y <- test.data$y
grp <- test.data$grp

weights <- rep(1/nrow(x), nrow(x))
sampleGrouping <- grp
covariateGrouping <- factor(1:ncol(x))
groupWeights <- c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping)))
parameterWeights <-  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x))
algorithm.config <- sgl.standard.config

# create data
data <- create.sgldata(x, y, weights, sampleGrouping)
sgl_test("sgl_test_dense", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, algorithm.config)

data <- create.sgldata(x, y, weights, sampleGrouping, sparseX = TRUE)
sgl_test("sgl_test_sparse", "sglOptim", data, covariateGrouping, groupWeights, parameterWeights, algorithm.config)
