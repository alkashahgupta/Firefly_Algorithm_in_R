#Reference Paper: Xin-She Yang, Nature-Inspired Metaheuristic Algorithms (second edition), pg 81-89, Luniver Press-2010

#You can either give your own input while running Firefly_algo in the bottom of this page or can leave default
#Inputs: no_firefly (total number of fireflies), N_iteration (number of timesteps),alpha (Randomness),betmin (Minimum value of beta (the attractiveness parameter)),gamma(Absorption coefficient), rng (minumum and maximum value), Fun (desired function to optimize)


firefly_algo <- function (no_firefly,N_iteration,alpha, betamin, gamma,rng,FunA = Fun){
	#Check input parameters (otherwise set as default values)
	# the input parameter is assumed in sequence
	nargin <- nargs()
	if (nargin < 1){
		rng <-c(-1,1)
		para <- c(20,50,0.25,0.20,1)
	}
	#Assume only one input parameter i.e. no_firefly is only passed to the function then rest of the parameter is initialized by the function
	else if (nargin<2){
		para <-c(no_firefly,50,0.25,0.20,1)
		rng <-c(-1,1)
	}
	else if (nargin < 3){
		para <-c(no_firefly,N_iteration,0.25,0.20,1)
		rng <-c(-1,1)
	}
	else if (nargin < 4){
		para <-c(no_firefly,N_iteration,alpha,0.20,1)
		rng <-c(-1,1)
	}
	else if (nargin < 5){
		para <-c(no_firefly,N_iteration,alpha, betamin,1)
		rng <-c(-1,1)
	}
	else if (nargin < 6){
		para <-c(no_firefly,N_iteration,alpha, betamin, gamma )
		rng <-c(-1,1)
	}
	else {
		rng<-c(rng[1],rng[2])
		para <-c(no_firefly,N_iteration,alpha, betamin,gamma)
	}
	
	#Assign values to each variables
	no_firefly <- para[1]
	MaxGeneration <- para[2]
	alpha <- para[3]
	betamin <- para[4]
	gamma <- para[5]
	
	
#Generate the initial locations of n fireflies
	
	init_ffa <- function (no_firefly,rng){
		ak<-rng[2] - rng[1]
		x1n <- t(matrix(runif(no_firefly)*ak + rng[1]))		
		x2n <- t(matrix(runif(no_firefly)*ak + rng[1]))
		Lightn <-matrix (0,nrow(x1n),ncol(x1n))
		return(list(x1n,x2n,Lightn))		#initial value before function evaluations
		
	}
 	
#Make sure the fireflies are within the range
 findrange_r <- function (x1n,x2n, rng){
 	for (i in 1: length(x1n)){
 		if (x1n[i] <= rng[1]){
 			x1n[i] <- rng[1]
 		}
 		if (x1n[i] >= rng[2]){
 			x1n[i] <- rng[2]
 		}
 		if (x2n[i]<= rng[1]){
 			x2n[i] <-rng[1]
 		}
 		if (x2n[i] >= rng[2]){
 			x2n[i]=rng[2]
 		}
 	}
 	return (list(x1n,x2n))
 }

		
	
#Move all fireflies to the better locations	
ffa_move <-function (x1n,x2n,Lightn,x1o,x2o,Lighto,alpha,gamma,rng){	
	ni = length(x1n)
	nj= length(x1o)
	for (i in 1:ni){
		for (j in 1:nj){
			r = sqrt((x1n[i]-x1o[j])^2 + (x2n[i]-x2o[j])^2) # This is the euclidian distance but you can change it according to ur needs
			
			if (Lightn[i]<Lighto[j]){
				beta0 =1
				#if you are not using betamin then beta <- beta0 * exp((-gamma) * r^2)
				beta<-(beta0-betamin)* exp(-gamma * r^2)+ betamin
				x1n[i] <- x1n[i]*(1-beta) + (x1o[j]*beta) +(alpha * (runif(1)-0.5)) # * exp(-0.2*log(downscale)))
				x2n[i] <- x2n[i]*(1-beta)+ (x2o[j] * beta) + (alpha * (runif(1) -0.5))	# * exp(-0.2*log(downscale)))
			}
		}
	}
	fdr <- findrange_r(x1n,x2n,rng);
	x1n <-fdr[[1]]
	x2n <-fdr[[2]]
	return (list(x1n,x2n))
	
}

		
## This function is optional, as it is not in the original FA. The idea to reduce randomness is to increase the convergence, however, if you reduce randomness too quickly, then premature convergence can occur. So use with care.
		
	alpha_new<-function(alpha,NGen){
		#alpha_n=alpha_0(1-delta)^NGen=0.005
		#alpha_0=0.9
		delta<-1-(0.005/0.9)^(1/NGen)
		alpha<-(1-delta)*alpha		
	}
	
#Initial values of an array
	zn <-rep((10^100),no_firefly)		
	As<-init_ffa(no_firefly,rng)
	x1n <-As[[1]]
	x2n <-As[[2]]
	Lightn <-As[[3]]

for (i in 1:MaxGeneration){	
	alpha=alpha_new(alpha,MaxGeneration) # This is optional
	prob <- c(0.8,0.2)
	zn <-FunA(x1n,x2n)
	# Ranking the fireflies by their light intensity
	tmp<-sort(zn,index.return=TRUE)
	Lightn<-tmp$x
	idx<-tmp$ix  #idx is index
	##Reordering of fireflies by their light intensity
	x1n<-x1n[idx] 
	x2n=x2n[idx]
	
	# initial value of xo, yo and Lighto is changed here	
	x1o<-x1n
	x2o<-x2n
	Lighto<-Lightn
	# Move all fireflies to the better locations
	tmp_move<- ffa_move(x1n,x2n,Lightn,x1o,x2o,Lighto,alpha,gamma,rng)
	x1n <-tmp_move[[1]]
	x2n <- tmp_move[[2]]
}
	# Output
	best <-cbind(x1n,x2n,Lightn)
	return (best)	
	
}

set.seed(1)
#Function to optimize (User might want to change this ) 
Fun <- function(x,y){
	str1=exp(-(x-4)^2-(y-4)^2)+exp(-(x+4)^2-(y-4)^2)+2*exp(-x^2-(y+4)^2)+2*exp(-x^2-y^2)
	return (str1)
}
firefly_algo()

