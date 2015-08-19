#read in output directory and library name
args=commandArgs(TRUE)

outdir=args[1]
lib=args[2]
sourcefile=args[3]
filename=args[4]

data=read.csv(sourcefile,header=FALSE,sep=" ")
tbl=t(cbind(data[2]))

jpeg(filename=paste(outdir,"/",lib,"_",filename,".jpg",sep='',collapse=''), width=800, height=800)

barplot(tbl, names.arg=t(data[1]), xlab="tag length (bp)", ylab="number of reads", border=NA, main=paste(lib, " - Position Along Index Trimmed Read With 3' Adapter"))

dev.off()




