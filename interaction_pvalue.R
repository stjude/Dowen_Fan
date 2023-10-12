options("scipen"=100, "digits"=4)

#=========================================================================================
# parameters
#=========================================================================================

args <- commandArgs()

print('THESE ARE THE ARGUMENTS')
print(args)

#ARGS
interaction_file = args[3]
peak_file = args[4]
FDR_cutoff= as.numeric(args[5])
PET_cutoff= as.numeric(args[6])
prefix = args[7]



#interaction_file="KZ358_filtered_PET_interactionSummary_succinct.txt"
#peak_file="KZ358_peaks_PETsummary_succinct.txt"
#FDR_cutoff=0.01
#PET_cutoff=3
#prefix="KZ358_filtered_PET"


#=========================================================================================
# calculation of p-value and FDR
#=========================================================================================


all_interactions=read.table(interaction_file,header=F,sep="\t",as.is=T,fill=T,comment.char = "#",quote="",skip=0,stringsAsFactor=F)
all_peaks=read.table(peak_file,header=F,sep="\t",as.is=T,fill=T,comment.char = "#",quote="",skip=0,stringsAsFactor=F)
all_interactions<-all_interactions[all_interactions$V1!=all_interactions$V2,]

cat("Number of interactions: ",nrow(all_interactions),"\n")
cat("Number of interactions with more than 1 PET: ",sum(all_interactions$V4>1),"\n")

index<-match(all_interactions$V1,all_peaks$V4)
interaction_5end<-all_peaks[index,]
index<-match(all_interactions$V2,all_peaks$V4)
interaction_3end<-all_peaks[index,]

## Calculate p-value
N=sum(all_interactions$V4)*2
p_values<-c()
p_values_byChr<-c()
for (i in 1:nrow(all_interactions)){
	if(i%%50000==0) print(i)
	i_ab=all_interactions$V4[i]	# number of PETs between two PET peaks
	c_a=interaction_5end$V5[i]	# number of PETs at 1st PET peak
	c_b=interaction_3end$V5[i]	# number of PETs at 2nd PET peak
	# p-value by hypergeometric distribution
	p_value=dhyper(i_ab, c_a, N-c_a, c_b)	
	p_values<-c(p_values,p_value)
	}

chr1<-sapply(all_interactions$V1, function(x) unlist(strsplit(x,":"))[1])
chr2<-sapply(all_interactions$V2, function(x) unlist(strsplit(x,":"))[1])
all_interactions$V6=ifelse(chr1==chr2,"intra","inter")
output<-data.frame(all_interactions,interaction_5end$V5,interaction_3end$V5,p_values)



## Calculate FDR for multiple hypothesis correction by 
bh_fdr=rep(-1,nrow(output))
bh_fdr<-p.adjust(p_values,method="fdr")
output$bh_fdr=bh_fdr

#length(bh_fdr)
#sum(output$bh_fdr<=0.01)
#sum(output$bh_fdr<=0.01 & output$V4>=3)


output_file<-paste(prefix,'_interactionSummary_bhFDR.txt',sep="")
colnames(output)<-c("clusterID_5end","clusterID_3end","interaction_ID","PETs","distance","class","clusterID_5end_count","clusterID_3end_count","p.value","BH.FDR")
write.table(output,file=output_file,row.names=F,col.names=T,sep="\t",quote=F) 


## intra-chromosomal interactions
output_file<-paste(prefix,'_interactionSummary_n',as.character(PET_cutoff),"FDR",as.character(FDR_cutoff),"_intrachromosomal.txt",sep="")
temp<-output[output$class=="intra" & output$BH.FDR<=FDR_cutoff & output$PETs>=PET_cutoff,]
write.table(temp,file=output_file,row.names=F,col.names=F,sep="\t",quote=F) 



chr1<-sapply(temp[,1], function(x) unlist(strsplit(x,":"))[1])
coor1<-sapply(temp[,1], function(x) unlist(strsplit(x,":"))[2])
coor2<-sapply(temp[,2], function(x) unlist(strsplit(x,":"))[2])
start_5end=sapply(coor1, function(x) as.numeric(unlist(strsplit(x,"-"))[1]))
end_5end=sapply(coor1, function(x) as.numeric(unlist(strsplit(x,"-"))[2]))
start_3end=sapply(coor2, function(x) as.numeric(unlist(strsplit(x,"-"))[1]))
end_3end=sapply(coor2, function(x) as.numeric(unlist(strsplit(x,"-"))[2]))

blockSizes=paste(as.character(end_5end-start_5end),as.character(end_3end-start_3end),sep=",")
blockStart=paste(rep("0,",nrow(temp)),as.character(start_3end-start_5end),sep="")

temp<-data.frame(chr1,start_5end,end_3end,temp[,3],temp[,4],"+",start_5end,end_3end,"0,0,0","2",blockSizes,blockStart)
output_file<-paste(prefix,'_interactionSummary_n',as.character(PET_cutoff),"FDR",as.character(FDR_cutoff),"_intrachromosomal.bed12",sep="")
write.table(temp,file=output_file,row.names=F,col.names=F,sep="\t",quote=F) 

## inter-chromosomal interactions
output_file<-paste(prefix,'_interactionSummary_n',as.character(PET_cutoff),"FDR",as.character(FDR_cutoff),"_interchromosomal.txt",sep="")
temp<-output[output$class=="inter" & output$BH.FDR<=FDR_cutoff & output$PETs>=PET_cutoff,]
write.table(temp,file=output_file,row.names=F,col.names=F,sep="\t",quote=F) 


## all high-confidence interactions in bedpe format
temp<-output[output$BH.FDR<=FDR_cutoff & output$PETs>=PET_cutoff,]
peak5end_index=match(temp$clusterID_5end,all_peaks$V4)
peak3end_index=match(temp$clusterID_3end,all_peaks$V4)
output_file<-paste(prefix,'_interactionSummary_n',as.character(PET_cutoff),"FDR",as.character(FDR_cutoff),".bedpe",sep="")
temp<-data.frame(all_peaks[peak5end_index,1:3],all_peaks[peak3end_index,1:3],temp[,3:4],".",".")
write.table(temp,file=output_file,row.names=F,col.names=F,sep="\t",quote=F) 


