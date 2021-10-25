
###  input parameters ####

# example vector of percentages
avector <- c(33, 28, 94, 73, 16, 71, 97, 98, 82, 94)
# m as group size (one gorup 2 mice, other 4)
m<-4
# The true element combi
true_combi<-c(1,2,3,4)

############## Get distribution ##########
n = length(avector)
combinations = combn(n, m)


subset_mean = function(x, vec){
  return(mean(vec[x]))}

subset_sum = function(x, vec){
  return(sum(vec[x]))}

subset_summed = apply(combinations, 2, subset_sum, vec = avector)
rest_mean <- (sum(avector)-subset_summed)/(n-m)
subset_mean = apply(combinations, 2, subset_mean, vec = avector)
diff_vector <- (subset_mean-rest_mean)


#### To get the p-value ##

true_value<-apply(combinations, 2, function(x) all(x == true_combi)) 
true_index<-min(which(true_value == TRUE))
true_diff<-diff_vector[true_index]

# Plot distribution with anything above true value coloured
dat <- data.frame( x=diff_vector, above=diff_vector>true_diff )
qplot(x,data=dat,geom="histogram",fill=above, binwidth=0.01)

## Calculate the p-value
m<-which(dat$above == FALSE)
P_value<- (1/n)*(1/(length(m)))
