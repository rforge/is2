library(is2)

pompExample(ou2)


p.truth <- coef(ou2)
guess2 <- guess1 <- p.truth
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.9*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.2*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
pdf("ou2mif2.pdf")
set.seed(64857673L)


is1a <- is2(ou2,Nmif=30,start=guess1,
             pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
             rw.sd=c(
               x1.0=.5,x2.0=.5,
               alpha.2=0.1,alpha.3=0.1),
             transform=F,
             Np=1000,
             var.factor=1,
             ic.lag=10,
             cooling.type="hyperbolic",
             cooling.fraction=0.05,
             method="ris1",
			 lag=1,
             tol=1e-8
             )
is2a <- is2(ou2,Nmif=30,start=guess1,
             pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
             rw.sd=c(
               x1.0=.5,x2.0=.5,
               alpha.2=0.1,alpha.3=0.1),
             transform=F,
             Np=1000,
             var.factor=1,
             ic.lag=10,
             cooling.type="hyperbolic",
             cooling.fraction=0.05,
             method="is2",
	     lag=1,
             tol=1e-8
             )


plot(c(is1a,is2a))

dev.off()
