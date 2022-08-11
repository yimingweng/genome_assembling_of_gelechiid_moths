
path <- getwd()

myhisto <- read.table(paste0(path, "/plin.histo"), header=F, sep="\t")
x11() # this only works for windows R I guess
plot(myhisto,
     type="l",
     xlab="kmer count",
     ylab="frequency") # plots the line from the data points 

# replot the histogram
highest <- myhisto[which(myhisto$V2 == max(myhisto[10:200,2])),]

x11() # this only works for windows R I guess
plot(myhisto[10:200,],
     type="l",
     xlab="kmer count",
     ylab="frequency") # plots the line from the data points 
points(myhisto[10:200,]) # plot the data points
text(highest[,1], highest[,2], highest$V1,
     cex=1, pos=2,col="black")