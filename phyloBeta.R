##a measure of phylobetadiversity (e.g., Graham & Fine 2008) using Lande's (1996) scheme, used phylocom format file (comID,abund,species name). Includes randomization test if n.rand > 0 (swap tip labels across the phylogeny).
#Marc Cadotte
#Nov 2008

PDcalc<- function (phy,com) {
	
	Community<-unique(com[,1])#generate list of community IDs
	N.coms<-length(Community)
	species<-phy$tip.label
	
	#for output
	PD<-numeric(N.coms)
	sp.rich<-numeric(N.coms)
	
	for (i in 1:N.coms) {
		
		ID<-Community[i]#select community.ID
		myclade<-com[com[,1]==ID,3]		
		if (length(myclade)==1) {
			print(paste("Community",ID,"has a single species!"))
		}
		
		sub.tree<-subTree(phy,myclade)
		
		PD[i]<-sum(sub.tree$edge.length)
		sp.rich[i]<-length(sub.tree$tip.label)
		
		
	}#end iter loop
	
	results <- data.frame(Community,PD,sp.rich)
	return(results)
}

subTree<- function(phy,species) {
	
	all.in<-is.na(match(species,phy$tip.label))
	
	if (sum(all.in)>0) 
		paste("species",as.character(species[all.in]),
			"not in phylogeny")
	
	dropme<-phy$tip.label[!phy$tip.label %in% unique(species)]
	sub.tr<-drop.tip(phy,dropme)
	
}

phyloBeta<-function (phy,com,n.rand=0) {
	
	PDT<-sum(phy$edge.length)
	
	PDa<-PDcalc(phy,com)
	
	PDa$q<-PDa$sp.rich/sum(PDa$sp.rich)	
	PDa$mean.PDa<-(PDa$q*PDa$PD)
	results<-data.frame(PD.alpha=sum(PDa$mean.PDa),
		PDT,PD.beta=PDT-sum(PDa$mean.PDa))
	
	if (n.rand>0) 	{
		
		tmp.phy<-phy
		rand.tmp<-data.frame(PD.alph.rand=0, PDT.rand=0,
		PD.beta.rand=0,n.rand=0)
		
		for (i in 1:n.rand) {
			
			tmp.phy$tip.label<-sample(tmp.phy$tip.label)
			rand.tmp[i,2]<-sum(tmp.phy$edge.length)
			tmp.PDa<-PDcalc(tmp.phy,com)
			tmp.PDa$q<-tmp.PDa$sp.rich/sum(tmp.PDa$sp.rich)
			rand.tmp[i,1]<-sum(tmp.PDa$q*tmp.PDa$PD)
			
			rand.tmp[i,3]<-rand.tmp$PDT.rand[i]-rand.tmp$PD.alph.rand[i]
			rand.tmp[i,4]<-n.rand	
			}
		results<-cbind(results,
		null.PD.alph=mean(rand.tmp$PD.alph.rand),
		SD.PD.alph=sd(rand.tmp$PD.alph.rand),
		null.PD.beta=mean(rand.tmp$PD.beta.rand),
		SD.PD.beta=sd(rand.tmp$PD.beta.rand),
		num.rands=n.rand)
		} 

	return(results)
	}