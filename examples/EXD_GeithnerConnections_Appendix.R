### Appendix to Geithner example
### J L'Hour
### 8 August 2018

### Gather old code that might be useful

### OTHER FIGURES

### B. Corrected
phi_qC = t(mapply(function(t) quantile(ResultP_C[,t], probs = c(.005,.025,.975,.995)), 1:length(TestPeriod)))
ATTdataC = ts(cbind(phi_qC[,1:2],phiP,phi_qC[,3:4]),start=c(-15), freq=1)

### Figure 3: Geithner connected firms effect vs. random permutations, corrected
pdf("plot/GeithnerAR_FisherTestCorrected.pdf", width=10, height=7)
plot(ATTdataC, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="AR, in pp",
     ylim=c(-.1,.1),
     main="")
abline(h=0,
       lty=2,col="grey")
lim <- par("usr")
rect(0, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(0.5,-.065,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()


ATTdata = ts(cbind(cumphi_q[,1:2],cumphiP,cumphi_q[,3:4]),start=c(1), freq=1)
### C. Figure 3: Geithner connected firms effect vs. random permutations
pdf("plot/GeithnerCAR_FisherTest.pdf", width=10, height=6)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="CAR, in pp",
     ylim=c(-.25,.25),
     main="Cumulative Abnormal Returns (CAR) for True Treatment vs. Random Permutations")
abline(h=0,
       lty=2,col="grey")
legend(1,-.15,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()

### CAR[0,1] and CAR[0,10] Table, corrected version

ATTdata = ts(cbind(cumphi_qC[,1:2],cumphiP,cumphi_qC[,3:4]),start=c(1), freq=1)
### D. Figure 4: Geithner connected firms effect vs. random permutations
pdf("plot/GeithnerCAR_FisherTestCorrected.pdf", width=10, height=6)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="CAR, in pp",
     ylim=c(-.25,.25),
     main="Cumulative Abnormal Returns (CAR) for True Treatment vs. Random Permutations, \n Corrected Inference")
abline(h=0,
       lty=2,col="grey")
legend(1,-.15,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()