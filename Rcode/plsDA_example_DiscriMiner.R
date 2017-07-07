library(DiscriMiner)
## Not run:
# load iris dataset
data(iris)

# PLS discriminant analysis specifying number of components = 2
my_pls1 = plsDA(iris[,1:4], iris$Species, autosel=FALSE, comps=2)
my_pls1$confusion
my_pls1$error_rate
# plot circle of correlations
plot(my_pls1)

# PLS discriminant analysis with automatic selection of components
my_pls2 = plsDA(iris[,1:4], iris$Species, autosel=TRUE)
my_pls2$confusion
my_pls2$error_rate

# linear discriminant analysis with learn-test validation
learning = c(1:40, 51:90, 101:140)
testing = c(41:50, 91:100, 141:150)
my_pls3 = plsDA(iris[,1:4], iris$Species, validation="learntest",
                learn=learning, test=testing)
my_pls3$confusion
my_pls3$error_rate
plot(my_pls3)
## End(Not run)