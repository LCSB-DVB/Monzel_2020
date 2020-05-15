###############################################################################
## (c) T.KAOMA, Luxembourg Institute of Health, 2019-08-27
## This is an adaptation of pvcaBatchAssess (pvca).
## The purpose of this is the usage of the function to bypass restriction
## for example, small number of feature ...
###############################################################################

suppressMessages({
    library("lme4")
    library(pvca)
    library('reshape')
})

pvca <- function (exprs, metadata, batch.factors, threshold, model.func){
	theDataMatrix = as.matrix(exprs)
	dataRowN <- nrow(theDataMatrix)
	dataColN <- ncol(theDataMatrix)
	theDataMatrixCentered = t(scale(t(exprs),center=T,scale=F))
	theDataCor <- cor(theDataMatrixCentered)
	eigenData  <- eigen(theDataCor)
	eigenValues = eigenData$values
	ev_n <- length(eigenValues)
	eigenVectorsMatrix = eigenData$vectors
	eigenValuesSum = sum(eigenValues)
	percents_PCs = eigenValues/eigenValuesSum
	# expInfo <- pData(abatch)[, batch.factors]
	expInfo <- metadata[,batch.factors]
	exp_design <- as.data.frame(expInfo)
	expDesignRowN <- nrow(exp_design)
	expDesignColN <- ncol(exp_design)
	# my_counter_2 = 0
	# my_sum_2 = 1
	# for (i in ev_n:1) {
		# my_sum_2 = my_sum_2 - percents_PCs[i]
		# if ((my_sum_2) <= threshold) {
			# my_counter_2 = my_counter_2 + 1
		# }
	# }
	my_counter_2 = which(cumsum(percents_PCs) >= threshold )[1]
	if (my_counter_2 < 3) {
		pc_n = 3
	}    else {
		pc_n = my_counter_2
	}
	# pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
		# pc_n), ncol = 1)
	# mycounter = 0
	# for (i in 1:pc_n) {
		# for (j in 1:expDesignRowN) {
			# mycounter <- mycounter + 1
			# pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
				# i]
		# }
	# }
	x = eigenVectorsMatrix[,1:pc_n]
	colnames(x) = paste("PC",1:pc_n,sep="")
	x = data.frame(ID=1:nrow(x),x)
	pc_data_matrix = matrix(melt(x,"ID")$value,ncol=1)
	rm(x)
	AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
	Data <- cbind(AAA, pc_data_matrix)
	variables <- c(colnames(exp_design))
	for (i in 1:length(variables)) {
		Data$variables[i] <- as.factor(Data$variables[i])
	}
	op <- options(warn = (-1))
    # effects_n = expDesignColN + choose(expDesignColN, 2) + 1
	## model.func <- c()
	## index <- 1
	## for (i in 1:length(variables)) {
		## mod = paste("(1|", variables[i], ")", sep = "")
		## model.func[index] = mod
		## index = index + 1
	## }
	#model.func = paste("(1|", variables, ")", sep = "")
	## for (i in 1:(length(variables) - 1)) {
		## for (j in (i + 1):length(variables)) {
			## mod = paste("(1|", variables[i], ":", variables[j], 
				## ")", sep = "")
			## model.func[index] = mod
			## index = index + 1
		## }
	## }
	# model.func = c(model.func,"(1|IDH_Mut:Treatment)","(1|IDH_Mut:Replicate)")
	effects_n = length(model.func) + 1
	function.mods <- paste(model.func, collapse = " + ")
	funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
	randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
	for (i in 1:pc_n) {
		y = (((i - 1) * expDesignRowN) + 1)
		Rm1ML <- lmer(
				funct, Data[y:(y-1+expDesignRowN), ],
				REML = TRUE, verbose = FALSE, 
				na.action = na.omit
			)
		randomEffects <- Rm1ML
		randomEffectsMatrix[i, ] <- c(
				unlist(VarCorr(Rm1ML)), 
				resid = sigma(Rm1ML)^2
			)
	}
	effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
    # randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
        # ncol = effects_n)
    # for (i in 1:pc_n) {
        # mySum = sum(randomEffectsMatrix[i, ])
        # for (j in 1:effects_n) {
            # randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                # j]/mySum
        # }
    # }
	randomEffectsMatrixStdze = randomEffectsMatrix/rowSums(randomEffectsMatrix)
	randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
        ncol = effects_n)
	for (i in 1:pc_n) {
		weight = eigenValues[i]/eigenValuesSum
		for (j in 1:effects_n) {
			randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
				j] * weight
		}
	}
	randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
	randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
	totalSum = sum(randomEffectsSums)
    randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
        ncol = effects_n)
    for (j in 1:effects_n) {
        randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
    }
	# return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
    pvcaObj <- list(dat = randomEffectsMatrixWtAveProp, label = effectsNames)
}
