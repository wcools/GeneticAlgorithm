# Genetic Algorithm
# for Design Of Experiments
# January, 2015
# w. cools

# genetic algorithms are used to search the parameterspace for an optimal solution
# various solutions (individuals) are looked into in parallel, each consisting of its building blocks (genes)
	# in current application the solution is an experiment, the building block is the condition under which the observation is made
# promising solutions are further explored others discarded based on a quality measure (fitness)
	# in current application the quality measure is the D-optimality (multivariate function that expresses the standard error of estimation)
	# note: D-optimality is determined for a model that includes variates of which the values are determined
# further exploration is done by recombining earlier solutions: recombination -> works with existing building blocks (genes)
# further exploration is done by altering earlier solutions: mutation -> creates new building blocks (genes)
	# in current application the current set of solutions is extended with an equally size set of recombined and mutated solutions
	# in current application the best (fittest) third of the solution is kept, the rest discarded

# note: genetic algorithms can be programmed in many different ways, some better than others and to some extent dependent on the purpose (no best algorithm exists)
	# note: in current application use is made of lists as datastructure which is conceptually clean but highly inefficient, use matrices instead


###############################
# < experiment = individual > #
###############################

# assumes fixed number of observations 
nrObs <- 16
# assume a set of parameters with boundaries
parBounds <- data.frame(t1=c(.3,1),t2=c(-1,1),t3=c(-1,1))
nrPars <- ncol(parBounds)
# a solution could be
     # t1   t2   t3
# 1  -2.8 -2.2  2.9
# 2  -2.8  2.5  1.1
# 3   0.6  3.0  0.3
# 4   1.0 -3.0  2.7
# 5  -2.7 -2.2  2.8
# 6   2.6  1.8 -2.8
# 7  -2.5  1.0  0.2
# 8  -1.0  2.9 -2.8
# 9   1.7 -2.2  1.9
# 10  3.0 -2.6 -2.5
# 11 -0.2 -1.7 -1.7
# 12 -2.8  2.9 -0.3
# 13 -2.6  3.0  2.8
# 14 -0.2  2.3  1.8
# 15  1.5 -2.8 -3.0
# 16  3.0 -1.8 -0.9

# create observations (genes) combined as experiment (individual)
setBoundObs <- function(nrObs){
	round(t(replicate(nrObs,runif(ncol(parBounds),as.numeric(parBounds[1,]),as.numeric(parBounds[2,])))),2)
}
checkConstraintObs <- function(mx){
	abs(rowSums(mx))>1 | abs(.5*mx[,1]+.3*mx[,2]+.2*mx[,3]) > 1 | unlist(lapply(1:3,function(.x) any(mx[,.x] < parBounds[1,.x]) | any(mx[,.x] > parBounds[2,.x])))
}
setConstraintObs <- function(mx){
	tmp <- mx[checkConstraintObs(mx),,drop=FALSE]
	while(nrow(tmp)>0){
	cat(nrow(tmp),"\n")
	mx[checkConstraintObs(mx),] <- setBoundObs(nrow(tmp))
	tmp <- mx[checkConstraintObs(mx),,drop=FALSE]
	}
	mx <- data.frame(mx)
	names(mx) <- names(parBounds)
	mx
}

#########################
# < quality = fitness > #
#########################

# assume a model
model <- ~1+t1+t2+t3+I(t1*t1)+t1:t2+t1:t3+t2:t3

# determine the quality of each solution: fitness
fitness <- function(modelMx){	return(det(t(modelMx)%*%modelMx))	}

#######################################
# < search/optimisation = evolution > #
#######################################
# a genetic algorithm searches the parameter-space for experiments
# new generations are created by generating new individuals

# the number of individuals within a population: parallel search
popSize <- 50
# the maximum number of generations: sequential search
nrGenerations <- 10

########################
# < first generation > #
########################

# a list of possible experiments is set up to start with, making sure that they are all in agreement with the constraints
pop <- vector("list",length=popSize)
pop <- lapply(pop, function(.x) .x <- setBoundObs(nrObs))
pop <- lapply(pop, function(.x) .x <- setConstraintObs(.x))

# evaluate all individuals within a population
evaluate <- function(pop){
	# keep parameters and add the modelled set of parameters
	out <- lapply(pop, function(.x) list(parameters=.x$parameters,model=model.matrix(model,.x$parameters)))
	# add the fitness and constraint check
	out <- lapply(out,function(.i) list(parameters=.i$parameters,model=.i$model,fit=fitness(.i$model),constraints=any(checkConstraintObs(.i$parameters))))
	# return a list of individuals (a generation)
	return(out)
}

# parents is list of individuals, each with following structure
parents <- evaluate(lapply(pop,function(.i) list(parameters=data.frame(.i))))
# List of 4
 # $ parameters :'data.frame':    16 obs. of  3 variables:
  # ..$ t1: num [1:16] 0.8 0.91 0.99 0.46 0.69 0.98 0.85 0.91 0.41 0.69 ...
  # ..$ t2: num [1:16] -0.26 0.14 -0.98 -0.09 0.1 -0.23 -0.57 -0.6 -0.06 -0.58 ...
  # ..$ t3: num [1:16] -0.58 -0.8 0.09 0.39 -0.78 -0.47 0.35 0.14 0.64 0.09 ...
 # $ model      : num [1:16, 1:8] 1 1 1 1 1 1 1 1 1 1 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:16] "1" "2" "3" "4" ...
  # .. ..$ : chr [1:8] "(Intercept)" "t1" "t2" "t3" ...
  # ..- attr(*, "assign")= int [1:8] 0 1 2 3 4 5 6 7
 # $ fit        : num 0.00442
 # $ constraints: logi FALSE

#######################
# < next generation > #
#######################

# parents (remain in pool) to avoid loosing proper solutions

# cross-over children
crossover <- function(parents,unit1,unit2){
	child <- parents[[unit1]]$parameters
	select <- sample(c(T,F),nrow(child),replace=T)
	child[select,] <- parents[[unit2]]$parameters[select,]
	return(list(parameters=child))
}

# mutate genes
mutate <- function(parents,unit,p,range){
	child <- parents[[unit]]$parameters
	Select <- matrix(sample(c(T,F),nrow(child)*ncol(child),prob=c(p,1-p),replace=T),ncol=ncol(child))
	tmp <- child
	while(sum(Select)>0){
		tmp[Select] <- round(child[Select]+runif(sum(Select),-1*range,range),2)
		if(!any(checkConstraintObs(tmp))){
			child <- tmp
			Select <- 0
		}
	}
	return(list(parameters=child))
}


#################
# < evolution > #
#################

# run over generations
g <- bestFit <- 0
nrGenerations <- 1000
probs <- seq(.05,.95,length.out=nrGenerations)
# for as long as...
while(g <= nrGenerations){

one <- Sys.time()

	combine <- matrix(sample(1:length(parents),prob=length(parents):1),ncol=2,byrow=T)
	children <- c(evaluate(apply(combine,1,function(.r) crossover(parents,.r[1],.r[2]))),
				evaluate(apply(combine,1,function(.r) crossover(parents,.r[2],.r[1]))))
	mutations <- evaluate(lapply(1:length(parents), function(.u) mutate(parents,.u,.1,probs[g+1])))
	population <- c(parents,children,mutations)
	population <- population[rev(order(unlist(lapply(population,function(.i) .i$fit))))]
	# population <- population[unlist(lapply(population,function(.i) .i$constraints))]
	parents <- population[1:length(parents)]

two <- Sys.time()

	g = g + 1
	cat(g,"> fit: ",parents[[1]]$fit," (",parents[[1]]$fit-bestFit,") taking: ",two-one,"\n")
	bestFit <- parents[[1]]$fit
	
	# bestFit <- c(bestFit,parents[[1]]$fit)
	# plot(bestFit)
	# if(var(unlist(lapply(population,function(.i) .i$fit)))[1:popSize]==0) g <- nrGenerations + 1
}


###############################
# < survival of the fittest > #
###############################

# final solution for the experiment
final <- parents[[1]]$parameters
     # t1    t2    t3
# 1  1.00 -0.48  0.48
# 2  0.67 -0.92 -0.75
# 3  1.00 -1.00  1.00
# 4  0.72  0.62 -0.34
# 5  1.00  0.99 -1.00
# 6  1.00  0.97 -1.00
# 7  1.00  0.99 -1.00
# 8  1.00 -1.00 -1.00
# 9  0.30  0.82 -0.12
# 10 0.65  0.93 -0.58
# 11 0.30  0.99 -1.00
# 12 0.30 -1.00  1.00
# 13 0.30  0.25  0.45
# 14 0.61 -1.00  1.00
# 15 0.30 -0.85 -0.45
# 16 0.30  1.00 -0.99
# note: while the best solution would be on the edges, it is not everywhere picked up as such, .99 should probably be 1

# for the resulting extended model
parents[[1]]$model
# with D-optimality value
parents[[1]]$fit

###############################

##### 3d plot of the best solution #####
library(scatterplot3d)
# for plotting: generate all possible points within the boundary
tmp <- apply(parBounds,2,function(.x) seq(.x[1],.x[2],length.out=30))
tmp <- expand.grid(tmp[,1],tmp[,2],tmp[,3])
tmp <- tmp[!(abs(rowSums(tmp))>1 | abs(.5*tmp[,1]+.3*tmp[,2]+.2*tmp[,3]) > 1),]
wer <- scatterplot3d(tmp[,1],tmp[,2],tmp[,3],type='n',xlab="t1",ylab="t2",zlab="t3",main=paste("best fit",as.character(model)))
wer$points3d(tmp[,1],tmp[,2],tmp[,3],col = "lightgray", type = "p", pch = 19,cex=.1)
wer$points3d(final,col = "red", type = "p", pch = 19,cex=-1*(final[,2]-1.5)/2+.5)
savePlot(file="BestFit3d.png",typ="png")

