source("http://www.bioconductor.org/biocLite.R")
biocLite("ALL")
library(ALL)
data(ALL)

#average gene expression per patient in 5 different visualizations
all.exprs <- exprs(ALL)
patient.means <- apply(all.exprs, 2, mean)
old.par <- par(mfrow = c(2, 3), ps = 20)
plot(hist(patient.means, plot = F))
boxplot(patient.means, main = "Average gene expression level per patient")
stripchart(patient.means, method = "jitter", main = "Average gene expression level per patient")
dotchart(patient.means, main = "Average gene expression level per patient")
plot(sort(patient.means), main = "Average gene expression level per patient")
par(old.par)

#The average gene expression values for each patient show relatively little variation across patients and resemble normally distributed values. 
#This is probably due to the fact that gene expression values were normalized to become comparable across patients. 

#effect of gene on age using ANOVA and running function across all 12000+ genes
anova.lm.pval <- function(x) {
    df.tmp <- data.frame(Expr = x, Age = pData(ALL)$age)
    anova(lm(Age ~ Expr, df.tmp))["Expr", "Pr(>F)"]
}
p.age <- apply(exprs(ALL), 1, anova.lm.pval)   
p.age[which.min(p.age)] #which gene has most significance and least significance based on p-value
p.age[which.max(p.age)]
#histogram of the p-values of all genes vs. age
hist(p.age, main = "Effect of gene expression on age", xlab = "Slope p-value",
       col = 'lightgreen')

#after finding the gene that had the most effect on aging, I graphed a more comprehensive relationship between age and gene expressions 
max.eff.gene.expr <- exprs(ALL)[which.min(p.age),]

plot(max.eff.gene.expr, pData(ALL)$age,
     xlab = "Gene Expression", ylab = "Age", main = names(which.min(p.age)))
min.max.expr <- c(min(max.eff.gene.expr), max(max.eff.gene.expr))
points(min.max.expr, predict(max.eff.gene.age.lm, data.frame(Expr = min.max.expr)),
      col = "red", type = "l", lwd = 2)

#created a linear model line to estimate gene expression based on age:
old.par <- par(mfrow = c(2, 2), ps = 20)
plot(max.eff.gene.age.lm)
par(old.par)



#similar to gene expression to age, there seemed to be a correlation between gene expression and time-to-remission
ALL.pdat <- pData(ALL)
date.cr.chr <- as.character(ALL.pdat$date.cr)
diag.chr <- as.character(ALL.pdat$diagnosis)
date.cr.t <- strptime(date.cr.chr, "%m/%d/%Y")   #manipulating data to accurately calculate length until cancer remission
diag.t <- strptime(diag.chr, "%m/%d/%Y")
days2remiss <- as.numeric(date.cr.t - diag.t)

p.value <- 0
for (i in 1:12625) {
    gene.expr <- exprs(ALL)[i,]
    my.df <- data.frame(gene = gene.expr, T2R = days2remiss)
    my.lm <- lm(T2R ~ gene, my.df)
    p.value[i] <- anova(my.lm)$'Pr(>F)'[1]
}

hist(p.value), xlab = 'P-value', main = 'Gene expression vs. Time to remission')

#Finding relation between gene expressions and cancer disease subtype (B or T cell cancer)

ALL.pdat <- pData(ALL)
bt <- factor(substring(ALL.pdat$BT, 1, 1)) # disease subtype represented as B or T 
# calculate p-value for each gene: 
bt.pvals <- apply(exprs(ALL), 1, function(x) anova(lm(x ~ bt))[1, "Pr(>F)"])
old.par <- par(mfrow = c(1, 2), ps = 20)
plot(hist(bt.pvals, plot = F), col = 'lightgreen') # plot the histogram of p-values 
# draw the boxplot of expression levels of the most significant gene, stratified # by disease subtype: 
boxplot(exprs(ALL)[which.min(bt.pvals),] ~ bt, main = rownames(exprs(ALL))[bt.pvals == min(bt.pvals)])
par(old.par)

#result showed that gene #38319_at had the most significant p-value for relation to subtype

