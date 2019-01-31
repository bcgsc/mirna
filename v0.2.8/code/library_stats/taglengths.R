#read in output directory and library name
args=commandArgs(TRUE)

outdir=args[1]
lib=args[2]
sourcefile=args[3]
filename=args[4]

data=read.csv(sourcefile,header=TRUE)

#combine data together to make graph clearer
other_rnas=data['snoRNA'] + data["tRNA"] + data["rRNA"] + data["snRNA"] + data["scRNA"] + data["srpRNA"] + data["rmsk_RNA"] + data["No_CDS"]
coding_gene=data['p5_UTR'] + data['p3_UTR'] + data["Coding_Exon"] + data["Intron"]
repeats=data['LINE'] + data['SINE'] + data['LTR'] + data['Satellite'] + data['rmsk_DNA'] + data['rmsk_Low_complexity'] + data['rmsk_Simple_repeat'] + data['rmsk_Other'] + data['rmsk_Unknown']


tbl=t(cbind(data["mature"], data["star"], data["unannotated"], data["crossmapped"], data["stemloop"], data["precursor"],
other_rnas,
coding_gene,
repeats,
data["Unknown"]))
colours=colors()[c(26, 563, 73, 435, 551, 598,
35,
38,
51,
325)]


jpeg(filename=paste(outdir,"/",lib,"_",filename,".jpg",sep='',collapse=''), width=800, height=800)

barplot(tbl, names.arg=t(data["taglen"]), col=colours, xlab="tag length (bp)", ylab="% of all reads aligning to <3 positions", border=NA, ylim=c(0,ceiling(max(colSums(tbl))/10)*10), main=paste(lib, " - Percentage of Aligned Tags At Each Tag Length With Annotation"))
legend("right",
c("miRNA, mature strand", "miRNA, star strand", "miRNA, unannotated in miRBase", "miRNA, crossmapped", "miRNA, stemloop", "miRNA, precursor",
"Other RNAs",
"Coding Genes",
"Repeats",
"Unannotated"),col=colours,lwd=3)

dev.off()




