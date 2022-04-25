#' Unmix sediment mixtures
#'
#' Asses the relative contribution of the potential sediment sources for each sediment mixture in the dataset.
#'
#' @param data Data frame containing sediment source and mixtures
#' @param samples Number of samples in each hypercube dimension
#' @param iter Iterations in the source variability analysis
#' @param Means Boolean to switch when using mean and sd data
#' @param seed Seed for the random number generator
#'
#' @return Data frame containing the relative contribution of the sediment sources for each sediment mixture and iterations
#'
#' @export
#'
unmix <- function(data, samples = 100L, iter = 100L, Means = F, seed = 123456L){
  system.time({
    if (Means == T) {
      sources <- data[c(1:(nrow(data)-1)),c(2:ncol(data))]
      x <- round((ncol(data)-3)/2 + 2)
      mixtures <- inputSample(data[nrow(data),1:x])
      colnames(sources)[1] <- "id"
      
    } else {
      ###########################################
      sources <- inputSource(data)
      mixtures <- inputSample(data)
      ###########################################
    }
    
    # verify the number of sources and properties
    if (nrow(sources) -1 >= ncol(mixtures) ) {
      warning("As a minimum, n - 1 properties are required to discriminate rigorously between n sources. Additional properties are frequently required to increase the reliability of the results")
    }
    
    nsources1<- as.data.frame(unique(data[,2]))
    nsources<- nsources1[-nrow(nsources1),] 
    nsources<- as.vector(nsources) 
    
    if (Means == T) {
      cat("Summary of the model imputs:
        ", (ncol(data)-3)/2, "variables from",nrow(nsources1)-1,"sources (",nsources,")",
          "\n")
      
    } else {
      cat("Summary of the model imputs:
        ", ncol(data)-2, "variables from",nrow(nsources1)-1,"sources (",nsources,")",
          "\n")
    }
    
    #invisible(readline(prompt="Press [enter] to unmix your data"))
    
    # results <- unmix_c(sources, mixtures, samples, 1, seed)  
    if (iter==1) {  
      results <- unmix_c(sources, mixtures, samples, iter, seed)
    }
    else {  
      results <- unmix_c(sources, mixtures, samples, iter+1, seed)
    }
    
    # reorder factor levels in order of appearance
    data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
    # read groups (second column)
    groups <- data[, 2]
    # asume last group is mixtures
    mixture <- levels(groups)[nlevels(groups)]
    # read sources
    sources <- data[!groups == mixture, ]
    # remove mixture level
    groups <- levels(droplevels(sources[, 2]))
    # replace column names
    colnames(results) <- c("id", "GOF", groups)
    
    # read groups (second column)
    groups <- data[, 2]
    # asume last group is mixtures
    mixture <- levels(groups)[nlevels(groups)]
    # read mixtures
    mixtures <- data[groups == mixture, ]
    # replace sample names
    results$id <- as.character(results$id)
    # count sources rows to modify mixture id in the results
    nrow_sources <- nrow(sources)
    
    for (i in 1:ncol(results)) {
      results[, i] <- as.numeric(as.character(results[, i]))
    }
    
    results$id <- results$id + nrow_sources
    
    
    results <- results[order(results[, 1]), ]
    rownames(results) <- 1:nrow(results)
    
    
    {
      if (iter==1) {  
        cat("Summary of the model outputs:",
            "\n",
            "See below the result/s of the unmixing process using the central value or the average with no correction",
            "\n",
            "\n")
        print(aggregate(. ~ id, data = results, function(x) c(mean = mean(x))))
      }  
      
      else {  
        cat("Summary of the model outputs:",
            "\n",
            "See below the result/s of the unmixing process using the source variability of the best", iter, "results, notice that the first row of the results is the central value or the average with no correction",
            "\n",
            "\n")
        print(aggregate(. ~ id, data = results, function(x) c(mean = mean(x), SD = sd(x))))
        
      }  
    }
    
    return(results)
    
  })
}

unmix_R_legacy <- function(source, mixture, iter = 2000, seed = 123456)
{
	set.seed(seed)

	snames <- source[,1]

	source <- data.matrix(source[-1])
	mixture <- data.matrix(mixture[-1])
	
	# normalize
	cols <- (ncol(source)-1)/2
	for (col in c(1:cols))
	{
		mx <- max(source[,col]+source[,cols+col])
		mn <- min(source[,col]-source[,cols+col])
		source[,col] <- ( source[,col] - mn ) / ( mx - mn )
		source[,cols+col] <- source[,cols+col] / ( mx - mn )
		mixture[,col] <- ( mixture[,col] - mn ) / ( mx - mn )
	}
	
	nsource <- nrow(source)
	ntracer <- (ncol(source)-1)/2

	r <- matrix(, nrow = iter, ncol = nsource+2)

	csource <- matrix(, nrow = nsource, ncol = ntracer)

	for (i in c(1:iter))
	{
		y <- c()
		x <- matrix(, nrow = ntracer, ncol = nsource-1)

		if(i==2 || i==3)
		{
			for (j in c(1:ntracer))
			{
				y <- c(y, mixture[1,j][[1]]-source[nsource,j])
				for (k in c(1:(nsource-1)))
				{
					x[j, k] <- source[k,j]-source[nsource,j]
				}
			}
			for (j in c(1:ntracer))
			{
				for (k in c(1:nsource))
				{
					csource[k,j] <- source[k,j]
				}
			}
		}
		else
		{
			ls <- c()
			for (j in c(1:ntracer))
			{
				x1 <- rt(1, source[nsource, ntracer*2+1]) / sqrt(source[nsource, ntracer*2+1])
				csource[nsource, j] <- source[nsource, j] + source[nsource, ntracer+j] * x1
				ls <- c(ls, csource[nsource, j])
			}
				
			for (j in c(1:ntracer))
			{
				y <- c(y, mixture[1,j][[1]] - ls[j])
				for (k in c(1:(nsource-1)))
				{
					x1 <- rt(1, source[k, ntracer*2+1]) / sqrt(source[k, ntracer*2+1])
					csource[k,j] <- source[k,j] + source[k,ntracer+j] * x1
					x[j, k] <- csource[k,j] - ls[j]
				}
			}
		}

		# least squares method
		model <- lm.fit(x, y)
		w <- as.vector(coef(model))
		w <- c(w, 1-sum(w))
		
		gof <- c()
		for (j in c(1:ntracer))
		{
			x1 <- 0.0
			for (k in c(1:nsource))
			{
				x1 <- x1 + w[k] * csource[k,j]
			}
			gof <- c(gof, (mixture[1,j][[1]]-x1)^2)
		}
		gof <- 1.0 - mean(gof)
		
		w <- c(1, gof, w)
		
		r[i, ] <- w
	}
	
	r <- as.data.frame(r)
	colnames(r) <- c('id', 'gof', paste( "w.", snames[c(1:nsource)], sep=""))
	return(r)
}


unmix_R <- function(source, mixture, iter = 2000, seed = 123456)
{
	set.seed(seed)

	snames <- source[,1]

	source <- data.matrix(source[-1])
	mixture <- data.matrix(mixture[-1])
	
	# normalize
	cols <- (ncol(source)-1)/2
	for (col in c(1:cols))
	{
		mx <- max(source[,col]+source[,cols+col])
		mn <- min(source[,col]-source[,cols+col])
		source[,col] <- ( source[,col] - mn ) / ( mx - mn )
		source[,cols+col] <- source[,cols+col] / ( mx - mn )
		mixture[,col] <- ( mixture[,col] - mn ) / ( mx - mn )
	}
	
	nsource <- nrow(source)
	ntracer <- (ncol(source)-1)/2

	r <- matrix(, nrow = iter, ncol = nsource+2)

	# compute central solution
	y <- c()
	x <- matrix(, nrow = ntracer, ncol = nsource-1)
	for (j in c(1:ntracer))
	{
		y <- c(y, mixture[1,j][[1]]-source[nsource,j])
		for (k in c(1:(nsource-1)))
		{
			x[j, k] <- source[k,j]-source[nsource,j]
		}
	}
	# least squares method
	model <- lm.fit(x, y)
	cw <- as.vector(coef(model))
	cw <- c(cw, 1-sum(cw))

	for (i in c(1:iter))
	{
		vm <- c(1:ntracer)*0		
		if(i==2 || i==3)
		{
			# use measured mixture
			for (j in c(1:ntracer))
			{
				vm[j] <- mixture[1,j][[1]]
			}
		}
		else
		{
			# compute virtual mixture
			for (j in c(1:ntracer))
			{
				vm[j] <- 0
				for (k in c(1:nsource))
				{
					x1 <- rt(1, source[k, ntracer*2+1]) / sqrt(source[k, ntracer*2+1])
					vm[j] <- vm[j] + cw[k] * ( source[k,j] + source[k,ntracer+j] * x1 )
				}
			}
		}

		y <- c()
		x <- matrix(, nrow = ntracer, ncol = nsource-1)
		for (j in c(1:ntracer))
		{
			y <- c(y, vm[j]-source[nsource,j])
			for (k in c(1:(nsource-1)))
			{
				x[j, k] <- source[k,j]-source[nsource,j]
			}
		}
		# least squares method
		model <- lm.fit(x, y)
		w <- as.vector(coef(model))
		w <- c(w, 1-sum(w))
		
		gof <- c()
		for (j in c(1:ntracer))
		{
			x1 <- 0.0
			for (k in c(1:nsource))
			{
				x1 <- x1 + w[k] * source[k,j]
			}
			gof <- c(gof, (vm[j]-x1)^2)
		}
		gof <- 1.0 - mean(gof)
		
		w <- c(1, gof, w)

		r[i, ] <- w
	}
	
	r <- as.data.frame(r)
	colnames(r) <- c('id', 'gof', paste( "w.", snames[c(1:nsource)], sep=""))
	return(r)
}

