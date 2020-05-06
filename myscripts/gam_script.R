require(mgcv)

gene6 <- dat_gam[,c(5, ncol(dat_gam))]

gam6 <- gam(V5 ~ s(treatment, k=6), family=nb(), data=gene6) # method doesn't seem to change edf so default is fine

plot(gam6, shade=TRUE, seWithMean=TRUE, scale=0, ylim=c(-4, 8))
points(x = conds, y=log2_q[9,], col='red')

gam.check(gam6)

predict(gam6)
gene6$V9

pd <- data.frame(treatment=seq(0, 3, .5))
plot(x=gene6$treatment, y=gene6$V5)
points(x=pd$treatment, y=exp(predict(gam6, newdata=pd)), col='red')
