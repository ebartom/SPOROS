

args <- commandArgs()

toxFile <-sub('--toxFile=', '', args[grep('--toxFile=', args)])

print("toxFile")
toxFile

toxes <- read.delim(file=toxFile,header=TRUE,sep="\t")
colnames(toxes) <- gsub(".txt.",".",colnames(toxes))
colnames(toxes) <- gsub(".geo.",".",colnames(toxes))

head(toxes)

png(paste(toxFile,"box.png",sep="."))
par(mai=c(3,0.82,0.4,0.42))
boxplot(toxes,las=2)
dev.off()
