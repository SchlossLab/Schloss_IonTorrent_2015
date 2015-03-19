################################################################################
#
# plot_error_profile.R
#
#
# Here we run the output of running seq.error through a script to generate a
# figure similar to what we presented as Figure 2 in Kozich et al. AEM 2013.
#
# Dependencies...
# * data/basic_enzyme?/mock?.trim.filter.error.matrix
# * data/basic_enzyme?/mock?.trim.filter.error.qual.reverse
# * data/basic_enzyme?/mock?.trim.filter.error.quality
# * data/basic_enzyme?/mock?.trim.filter.error.seq.reverse
# * data/basic_enzyme?/mock?.trim.filter.error.summary
# * data/basic_enzyme?/mock?.trim.summary
#
# Produces...
# * results/figures/basic_error_profile.png
#
################################################################################


figure <- "results/figures/basic_error_profile.png"
stub <- paste0("data/basic_enzyme")

#png(file=figure)

par(mfcol=c(2,2))

#A
par(mar=c(1, 5, 5, 1))
a <- read.table(file="data/basic_enzyme1/mock1.trim.filter.error.seq.reverse", header=T)
b <- read.table(file="data/basic_enzyme1/mock2.trim.filter.error.seq.reverse", header=T)
c <- read.table(file="data/basic_enzyme1/mock3.trim.filter.error.seq.reverse", header=T)

plot(100*(1-a$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)
plot(100*(1-b$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)
plot(100*(1-c$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)

nr <- min(c(nrow(a), nrow(b), nrow(c)))
composite <- (a[1:nr,] + b[1:nr,] + c[1:nr,])/3

#plot(100*(1-composite$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)
plot(100*(1-c$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)

a <- read.table(file="data/basic_enzyme2/mock1.trim.filter.error.seq.reverse", header=T)
b <- read.table(file="data/basic_enzyme2/mock2.trim.filter.error.seq.reverse", header=T)
c <- read.table(file="data/basic_enzyme2/mock3.trim.filter.error.seq.reverse", header=T)

plot(100*(1-a$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)
plot(100*(1-b$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)
plot(100*(1-c$match[1:nr]), xlim=c(0,nr), ylim=c(0,30), type="l", xlab="", ylab="", xaxt="n", cex.lab=1.2)



#Substitution rate (%)





#C
par(mar=c(1, 5, 5, 1))
a <- read.table(file="data/basic_enzyme1/mock1.trim.filter.error.quality"), header=T, row.names=1)
b <- read.table(file="data/basic_enzyme1/mock2.trim.filter.error.quality"), header=T, row.names=1)
c <- read.table(file="data/basic_enzyme1/mock3.trim.filter.error.quality"), header=T, row.names=1)
composite <- a + b + c

q <- as.numeric(rownames(composite))

m<-sample(rep(q, composite$matches), 10000)
boxplot(m, at=1, xlim=c(0.5,4.5), ylim=c(0,40), ylab="Quality score", cex.lab=1.2)
boxplot(rep(q, composite$substitutions), at=2,add=T, outline=F, yaxt="n")
boxplot(rep(q, composite$insertions), at=3, add=T, outline=F, yaxt="n")
boxplot(rep(q, composite$ambiguous), at=4, add=T, outline=F, yaxt="n")
text(0.5, 39.9, label="C", cex=1.5, font=2)


#D
par(mar=c(6, 5, 0, 1))
a<-read.table(file=paste0(stub, "/Mock1_S1_L001_R2_001.rc.filter.error.quality"), header=T)
b<-read.table(file=paste0(stub, "/Mock2_S2_L001_R2_001.rc.filter.error.quality"), header=T)
c<-read.table(file=paste0(stub, "/Mock3_S3_L001_R2_001.rc.filter.error.quality"), header=T)
composite <- a + b + c

q <- as.numeric(rownames(composite))

m<-sample(rep(q, composite$matches), 10000)
boxplot(m, at=1, xlim=c(0.5,4.5), ylim=c(0,40), ylab="")
boxplot(rep(q, composite$substitutions), at=2,add=T, outline=F, yaxt="n")
boxplot(rep(q, composite$insertions), at=3, add=T, outline=F, yaxt="n")
boxplot(rep(q, composite$ambiguous), at=4, add=T, outline=F, yaxt="n")
axis(1, at=c(1,2,3,4), labels=c("Matches", "Substitutions", "Insertions", "Ambiguous"), las=2)
text(0.5, 39.9, label="D", cex=1.5, font=2)
par(mfrow=c(1,1))

dev.off()
