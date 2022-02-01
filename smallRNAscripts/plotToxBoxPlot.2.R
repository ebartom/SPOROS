library(ggplot2)

args <- commandArgs()

toxFile <-sub('--toxFile=', '', args[grep('--toxFile=', args)])

print("toxFile")
toxFile

toxes <- read.delim(file=toxFile,header=TRUE,sep="\t")
colnames(toxes) <- gsub(".txt.",".",colnames(toxes))
colnames(toxes) <- gsub(".geo.",".",colnames(toxes))

head(toxes)

numDown <- length(na.omit(toxes$delta.dn.txt))
numUp <- length(na.omit(toxes$delta.up.txt))

print(numDown)
print(numUp)

toxes2 <- data.frame(toxvalues = c(toxes$delta.dn.txt[1:numDown],toxes$delta.up.txt[1:numUp]),
                     toxtype =c("delta DN", "delta UP"))

toxes2$toxtype[1:numDown] <- "delta DN"
toxes2$toxtype[((numDown+1):(numDown+numUp))] <- "delta UP"

png(paste(toxFile,"box.png",sep="."))
par(mai=c(3,0.82,0.4,0.42))
ggplot(toxes2, aes(y=toxvalues, x=toxtype)) +
	geom_boxplot() +
  	stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, linetype = "solid", col = "blue", size=1.2) +
  	stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, linetype = "solid", col = "red", size=1.2 ) +
	annotate(x=2.4,y=3,geom="text",label="mean", col = "red" ) +
	annotate(x=2.4,y=6,geom="text",label="median", col = "blue" )	
		 
dev.off()
