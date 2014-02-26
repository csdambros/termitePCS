library("ape")

isoptera<-read.csv("Isoptera.csv",row.names=1)




sp.ls<-unique(isoptera[,c("Taxon","Genus","Subfamily","Family","Order")])

# Testa se uma mesma espécie está em mais de um gênero, subfamília, etc
sp.ls[!is.na(match(sp.ls[,1],sp.ls[duplicated(sp.ls[,1]),1])),][order(sp.ls[!is.na(match(sp.ls[,1],sp.ls[duplicated(sp.ls[,1]),1])),][,1]),]

###

sp.ls.int<-as.matrix(data.frame(lapply(sp.ls,as.numeric)))

#sp.ls.int[,1]<-0

#rownames(sp.ls.int)<-sp.ls[,1]

sp.ls.int[is.na(sp.ls.int)]<-0

#ver<-t(t(sp.ls.int)/apply(sp.ls.int,2,max,na.rm=TRUE))

sp.ls.dist<-(sp.ls.int[,1])%*%t(sp.ls.int[,1]*0+1)!=(sp.ls.int[,1]*0+1)%*%t(sp.ls.int[,1])

for(i in 2:ncol(sp.ls.int)){
	sp.ls.dist<-
		as.matrix(sp.ls.dist)+
		as.matrix((sp.ls.int[,i])%*%t(sp.ls.int[,i]*0+1)!=(sp.ls.int[,i]*0+1)%*%t(sp.ls.int[,i]))
}

dimnames(sp.ls.dist)<-list(sp.ls[,1],sp.ls[,1])

write.table(sp.ls.dist,"taxodist.csv")

termite.taxophylo<-as.phylo(hclust(as.dist(sp.ls.dist)))

termite.taxophylo$tip.label<-1:232

write.tree(termite.taxophylo,"termite.taxophylo.tre")
