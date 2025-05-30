setwd("")

#load your packages. you first need to install them.
library(phytools)
library(ape)
library(nlme)
library(geiger)
library(ggplot2)
library(devtools)
library(reshape)
library(calibrate)
library(qpcR)


##Testing the ratio data from the chela (ratio of Zn/Carbon at average of major & minor tooth) 
#and telson for Zn/Fe/Mn/Ca (ratio of Zn/Carbon at position 2,3,4)


#Read the Santibanez et al 2022 independend dated phylogeny
tree<-read.tree("santibanez2022_mcmc_independent_timetree.tre")
plot(tree)


#Read in your data on metal ions
metal.ions_check = read.csv("Test_1.csv", row.names = 1)
head(metal.ions_check)
metal.ions = read.csv("Test_1_8.csv", row.names = 1)
head(metal.ions)

#check if data matches tree
obj <- name.check(tree,metal.ions_check)
obj

#or check if species in the tree are in the database
setdiff(tree$tip.label, row.names(metal.ions))

# Prune out the species that are not in your data
trimmed.tree <- drop.tip(tree,obj$tree_not_data)
plot(trimmed.tree)

#rename the labels of the tree to fit
print(trimmed.tree$tip.label)
trimmed.tree$tip.label <- c("lvar" = "Lychas mucronatus","Htrill"="Hottentotta hottentotta","Bisra"="Buthus occitanus","Aaustralis"="Androctonus australis",
                            "LquiE"="Leiurus quinquestriatus","Parabuthus"="Parabuthus granulatus","Uvit" = "Uroplectes vittatus","Tarc"="Tityus trinitatis",
                            "Ccar"="Centruroides suffusus","Bburmeiste"="Bothriurus bonariensis","Brotheas"= "Brotheas sp","Hadr"= "Hadrurus arizonesis", "Pbae"="Paruroctonus boreus",
                            "Pspinig" ="Paravaejovis spinigerus", "Ddiablo"="Diplocentrus lindo","Pandinus"= "Pandinus imperator","Hpaucidens"= "Hadogenes phylladus","Opisasper"= "Opisthcanthus lepturus")

print(trimmed.tree$tip.label)
plot(trimmed.tree)
# now check that everything matches
name.check(trimmed.tree,metal.ions)

#sort our trait data to match the order of the tips in the tree. Many phylogenetic comparative method functions in R assume that data are arranged in an order that matches the storage of the tips of the phylogenetic tree at hand, so I always do this just in case.
	
	metal.ions[1,1] #in first row, first column
	metal.ions[1,] #first row
	metal.ions<-metal.ions[trimmed.tree$tip.label,] #tells it you want it in this specific order
	
	head(metal.ions)
	
# turn your data into vectors
	
	#P_zinc = major tooth Zinc conc (normalised) vs P_zinc_raw = not normalised data
	#T_Zinc = Position 2 telson Zinc conc
	#P_iron = Major tooth Iron conc
	
	P_zinc <- metal.ions[,1]	# Make a vector of Pedipalp Zinc
	names(P_zinc) <- row.names(metal.ions) # Name the vector elements.
	head(P_zinc)
	
	P_iron <- metal.ions[,2]
	names(P_iron) <- row.names(metal.ions)
	head(P_iron)
	
  P_calcium <- metal.ions[,3]
	names(P_calcium) <- row.names(metal.ions)
	head(P_calcium)
	
	T_zinc <- metal.ions[,4]	# Make a vector of Telson Zinc
	names(T_zinc) <- row.names(metal.ions)
	head(T_zinc)
	
	T_Mn <- metal.ions[,5]	# Make a vector of Telson Manganese
	names(T_Mn) <- row.names(metal.ions)
	head(T_Mn)
	
	T_calcium <- metal.ions[,6]	# Make a vector of Telson Calcium
	names(T_calcium) <- row.names(metal.ions) 
	head(T_calcium)
	
	P_CAR <- metal.ions[,8]  #Make a vector of Pedipalp Chela aspect ratio (CAR) 
	#7 is old ratio, 8 is new log10 ratio
	names(P_CAR) <- row.names(metal.ions)
	head(P_CAR)
	
	###continue the same for the individual measurements if you want . . . 
	#etc.
	
	#ancestral reconstructon for P_zinc
	P_zinc.anc <- fastAnc(trimmed.tree,P_zinc,CI=TRUE)
	P_zinc.anc
	
	P_iron.anc <- fastAnc(trimmed.tree,P_iron,CI=TRUE)
	P_iron.anc
	
	P_calcium.anc <- fastAnc(trimmed.tree,P_calcium,CI=TRUE)
	P_calcium.anc
	
	T_zinc.anc <- fastAnc(trimmed.tree,T_zinc,CI=TRUE)
	T_zinc.anc
	
	T_Mn.anc <- fastAnc(trimmed.tree,T_Mn,CI=TRUE)
	T_Mn.anc
	
	T_calcium.anc <- fastAnc(trimmed.tree,T_calcium,CI=TRUE)
	T_calcium.anc

	P_CAR.anc <- fastAnc(trimmed.tree, P_CAR, CI=TRUE)
	P_CAR.anc

	
	#plot pedipalp and telson for zinc, side by side
	parsettings <- par(no.readonly=TRUE) # Record plotting parameters
	par(mfrow=c(1,2), cex=0.4)
	
	Chela_Zn <- contMap(trimmed.tree,P_zinc,res=30,fsize=c(2,0.5))
	Chela_Fe <- contMap(trimmed.tree,P_iron,res=30,fsize=c(2,0.5))
	Chela_Ca <- contMap(trimmed.tree,P_calcium,res=30,fsize=c(2,0.5))
	Chela_CAR <- contMap(trimmed.tree, P_CAR, res=30, fsize=c(2,0.5))
	Telson_Zn <- contMap(trimmed.tree,T_zinc,res=30,fsize=c(2,0.5))
	Telson_Mn <- contMap(trimmed.tree,T_Mn,res=30,fsize=c(2,0.5))
	Telson_Ca <- contMap(trimmed.tree,T_calcium,res=30,fsize=c(2,0.5))
	
		par(parsettings) # Reset to original plotting settings
	

	
	########test which model of evolution is appropriate
	#for P_zinc vs T_zinc
	# dependent variable is P_zinc and independent variable is T_zinc
	
	##### PEDIPALP (Zn) VS TELSON (Zn) COMPARISON#####
		
	result.lambda<-gls(P_zinc~T_zinc, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Zn~metal.ions$Telson_R_Zn)
	abline(result.lambda)
	
	result.ACDC<-gls(P_zinc~T_zinc, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Chela_R_Zn~metal.ions$Telson_R_Zn)
	abline(result.ACDC)
	
	result.brownian<-gls(P_zinc~T_zinc, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Chela_R_Zn~metal.ions$Telson_R_Zn)
	abline(result.brownian)
	
	result.OU<-gls(P_zinc~T_zinc, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Chela_R_Zn~metal.ions$Telson_R_Zn)
	abline(result.OU)
	
	##### TELSON (Ca) vs PEDIPALP (Ca) ####
	
	result.lambda<-gls(T_calcium~P_calcium, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
  plot(metal.ions$Telson_R_Ca~metal.ions$Chela_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(T_calcium~P_calcium, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.ACDC)
	
	result.brownian<-gls(T_calcium~P_calcium, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.brownian)
	
	result.OU<-gls(T_calcium~P_calcium, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.OU)
	
	##### PEDIPALP (Ca) VS PEDIPALP (Fe) COMPARISON #####
	
	result.lambda<-gls(P_calcium~P_iron, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Ca~metal.ions$Telson_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(P_calcium~~P_iron, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.ACDC)
	
	result.brownian<-gls(P_calcium~~P_iron, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.brownian)
	
	result.OU<-gls(P_calcium~P_iron, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.OU)
	
	
	##### PEDIPALP (Ca) VS PEDIPALP (Zn) COMPARISON #####
	
	result.lambda<-gls(P_calcium~P_zinc, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Ca~metal.ions$Telson_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(P_calcium~P_zinc, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.ACDC)
	
	result.brownian<-gls(P_calcium~P_zinc, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.brownian)
	
	result.OU<-gls(P_calcium~P_zinc, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Chela_R_~metal.ions$Telson_R_)
	abline(result.OU)
	
	
	
	##### TELSON (Zn) vs TELSON (Mn)COMPARISON ####
	
	result.lambda<-gls(T_zinc~T_Mn, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.lambda)
	
	result.ACDC<-gls(T_zinc~T_Mn, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.ACDC)
	
	result.brownian<-gls(T_zinc~T_Mn, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.brownian)
	
	result.OU<-gls(T_zinc~T_Mn, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.OU)
	
	##### TELSON (Zn) VS TELSON (Ca) COMPARISON ####
	
	result.lambda<-gls(T_zinc~T_calcium, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(T_zinc~T_calcium, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.ACDC)
	
	result.brownian<-gls(T_zinc~T_calcium, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.brownian)
	
	result.OU<-gls(T_zinc~T_calcium, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.OU)
	
	
	##### TELSON (Mn) VS TELSON (Ca) COMPARISON ####
	
	
	result.lambda<-gls(T_Mn~T_calcium, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Telson_R_Mn~metal.ions$Telson_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(T_Mn~T_calcium, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.ACDC)
	
	result.brownian<-gls(T_Mn~T_calcium, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.brownian)
	
	result.OU<-gls(T_Mn~T_calcium, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn_N~metal.ions$Telson_R_Mn_N)
	abline(result.OU)
	
	
	
	
	##### PEDIPALP (Zn) VS PEDIPALP (Fe) COMPARISON#####
	
    	result.lambda<-gls(P_zinc~P_iron, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Zn~metal.ions$Chela_R_Fe)
	abline(result.lambda)
	
	result.ACDC<-gls(P_zinc~P_iron, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.ACDC)
	
	result.brownian<-gls(P_zinc~P_iron, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.brownian)
	
	result.OU<-gls(P_zinc~P_iron, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.OU)
	##### PEDIPALP (Fe) VS PEDIPALP (Ca) COMPARISON ####
	result.lambda<-gls(P_iron~P_calcium, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Fe~metal.ions$Chela_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(P_iron~P_calcium, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.ACDC)
	
	result.brownian<-gls(P_iron~P_calcium, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.brownian)
	
	result.OU<-gls(P_iron~P_calcium, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.OU)
	
	##### PEDIPALP (Zn) VS PEDIPALP (Ca) COMPARISON ####
	result.lambda<-gls(P_zinc~P_calcium, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Zn~metal.ions$Chela_R_Ca)
	abline(result.lambda)
	
	result.ACDC<-gls(P_zinc~P_calcium, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.ACDC)
	
	result.brownian<-gls(P_zinc~P_calcium, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.brownian)
	
	result.OU<-gls(P_zinc~P_calcium, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.OU)
	
	##### PEDIPALP (Zn) VS PEDIPALP (CAR) COMPARISON ####
	result.lambda<-gls(P_zinc~P_CAR, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Zn~metal.ions$CAR_ratio)
	abline(result.lambda)
	
	result.ACDC<-gls(P_zinc~P_CAR, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.ACDC)
	
	result.brownian<-gls(P_zinc~P_CAR, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.brownian)
	
	result.OU<-gls(P_zinc~P_CAR, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.OU)
	
	
	
	##### PEDIPALP (Fe) VS PEDIPALP (CAR) COMPARISON ####
	result.lambda<-gls(P_iron~P_CAR, data=metal.ions, correlation=corPagel(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.lambda)
	plot(metal.ions$Chela_R_Fe~metal.ions$CAR_ratio)
	abline(result.lambda)
	
	result.ACDC<-gls(P_iron~P_CAR, data=metal.ions, correlation=corBlomberg(value=1, phy=trimmed.tree, fixed=TRUE),method="ML")
	summary(result.ACDC)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.ACDC)
	
	result.brownian<-gls(P_iron~P_CAR, data=metal.ions, correlation=corBrownian(value=1, phy=trimmed.tree),method="ML")
	summary(result.brownian)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.brownian)
	
	result.OU<-gls(P_iron~P_CAR, data=metal.ions, correlation=corMartins(value=1, phy=trimmed.tree, fixed=FALSE),method="ML")
	summary(result.OU)
	#plot(metal.ions$Telson_R_Zn~metal.ions$Telson_R_Mn)
	abline(result.OU)
	
	
	######TEST FOR BEST FIT REGRESSION MODEL#####
	
	
	regression.models_zinc<-matrix(,4,5,dimnames = list(c("Brownian Motion", "Early Burst", "Ornstein-Uhlenbeck", "Lambda"),c("log likelihood", "AIC", "Delta AIC", "AIC Weights", "Slope")))
	
	regression.models_zinc[,1]<-c(logLik(result.brownian)[[1]], logLik(result.ACDC)[[1]], logLik(result.OU)[[1]], logLik(result.lambda)[[1]])
	regression.models_zinc[,2]<-c(AIC(result.brownian), AIC(result.ACDC), AIC(result.OU), AIC(result.lambda))
	aic.all_zinc<-as.matrix(regression.models_zinc[,2])	
	scor.wts_zinc<-akaike.weights(aic.all_zinc)
	regression.models_zinc[,3]<-scor.wts_zinc$deltaAIC
	regression.models_zinc[,4]<-scor.wts_zinc$weights
	regression.models_zinc[,5]<-c(result.brownian[[4]][2], result.ACDC[[4]][2], result.OU[[4]][2], result.lambda[[4]][2])
	
	regression.models_zinc
	
	#sort by AIC score
	regression.models_zinc[order(regression.models_zinc[,2]),]
	

