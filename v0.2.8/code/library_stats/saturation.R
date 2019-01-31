#read in output directory and library name
args=commandArgs(TRUE)

outfile=args[1]
name=args[2]
sourcefile=args[3]

data=read.csv(sourcefile,header=TRUE)
data['Total.miRNA'] = data['Total.miRNA'] / 1000000
sat1x=cbind(data['Total.miRNA'],data['miRNA.Species'])
#sat10x=cbind(data['Total.miRNA'],data['miRNA.Species.Covered.by....10.Reads'])
sat10x=cbind(data['Total.miRNA'],data['miRNA.Species.Covered.by....10.Reads'])

jpeg(filename=outfile, width=800, height=800)
par(cex = 1.5)
plot(sat1x, col="blue", xlab="# reads aligned to miRNAs (Millions)", ylab="# miRNA species", ylim=c(0, ceiling(max(sat1x['miRNA.Species'])/100)*100), xlim=c(0, (ceiling(max(sat1x['Total.miRNA'])*1000000/250000)*250000) / 1000000), main=paste("miRNA Saturation in", name))
points(sat10x, col="red")
legend("bottomright", c(">= 1x coverage", ">= 10x coverage"), col=c("blue", "red"), pch=1, cex=.75)

#draw line of best fit if there are many samples
if (nrow(sat1x) > 4) {
	lines(smooth.spline(sat1x,df=4),col="blue",lwd=1)
	lines(smooth.spline(sat10x,df=4),col="red",lwd=1)
}

dev.off()
