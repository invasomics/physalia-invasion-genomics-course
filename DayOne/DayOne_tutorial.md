# Day One

Today, we'll be working with the Qfly dataset. We'll use four key methods to examine population structure and introduction pathways:
1. Population differentiation metrics (FST)
2. Population structure (PCA)
3. Population admixture (LEA)
4. Population assignment (DAPC)

Let's get into it!

### Open R and load in required packages
First, set your working directory to something logical (e.g., a 'DayOne' folder on your desktop), and run all code from within that directory.
```
library(vcfR)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(corrplot)
library(reshape2)
```
Depending on the version of R you are using, you may need to install some of these packages first. Do this with the following:
```
install.packages("vcfR")
install.packages("adegenet")
install.packages("StAMPP")
```
Then re-run the relevant library command, for example:
```
library(vcfR)
```
### Read in the Qfly dataset
Download the Qfly.vcf and Qfly_popmap.txt files from github (DayOne/RequiredFiles/) and put them into your working directory, then read them into R:
```
snp_vcf2 <- read.vcfR("Qfly.vcf")
pop.data2 <- read.table("Qfly_popmap.txt", header = F)
gl.snp2 <- vcfR2genlight(snp_vcf2)
pop(gl.snp2) <- rep(pop.data2$V2)
```
A few things to note:
1. The dataset is now a 'Genlight' object, with 6,707 SNPs and 310 individuals
2. There are a total of 28 populations in column 2 (V2) of the popmap file; each of the population names includes a prefix of 'Na', 'Ex', or 'In' (e.g., In1, In2) to show their native, expanded, and invasive status. I have also included a third column (V3), where I've reduced the dataset to three populations, corresponding to 'Native', 'Expanded', and 'Invasive'. You can confirm this by typing:
```
pop.data2$V2
```
and
```
pop.data2$V3
```

### FST
Let's calculate pairwise FST for all populations using the StAMPP package (the first line will take a few minutes to run). You can see more on the stamppFst command here: https://www.rdocumentation.org/packages/StAMPP/versions/1.6.3/topics/stamppFst (note nclusters = 6 sets the number of proccesor treads or cores to use during calculations and may need to be adjusted for your system):
```
Qfly_Fst <- stamppFst(gl.snp2, nboots = 100, percent = 95, nclusters = 6)
Fst <- Qfly_Fst$Fsts
pFst <- Qfly_Fst$Pvalues
write.table(Fst, "Fst.txt", sep="\t")
write.table(pFst, "Fst_pvalue.txt", sep="\t")
```
In this command, we are applying the stamppFst function to the gl.snp2 genlight object. This code produces two output files - one is a matrix of the pairwise Fst values and the other is a matrix of Fst p-values (i.e., telling us whether the differentiation is statistically significant or not). 

Let's create a heatmap of the FST values using the 'melt' function from R's reshape2 package:
```
Q_Fst <- as.matrix(read.table("Fst.txt"))
QflyFs <- melt(Q_Fst, na.rm = TRUE)
summary(QflyFs$value)
```
And plot:
```
ggplot(data = QflyFs, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#ffd60a", mid = "#4e9de6", high = "#001d3d", 
                       midpoint = 0.056, limit = c(0.005,0.11), space = "Lab", 
                       name="Pairwise Fst") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
```
In this plot, I've used a 3-colour gradient for the FST values and the populations are labelled.
Examine the plot and describe its key trends with respect to the three main populations we have (i.e., native, expanded, invasive). Take a moment to think about what you might expect this plot to look like and whether it meets your expectations.

You should see a pattern of higher differentiation among the pairwise comparisons that include the invasive populations, and lower values for the comparisons involving expanded and native populations. This is consistent with the admixture results (see below).

### PCA
Now let's use the adgenet package to perform a PCA analysis (see more on the glPca command here: https://www.rdocumentation.org/packages/adegenet/versions/2.1.11/topics/glPca):
```
snp.pca2 <- glPca(gl.snp2, nf = 10)
```
Write PCA scores to a file:
```
snp.pca.scores2 <- as.data.frame(snp.pca2$scores)
snp.pca.scores2$pop <- pop(gl.snp2)
write.table(snp.pca.scores2, "Qfly_adegenetPCA.txt", sep = "\t")
```
Export a list of eigen values and percentage variances for each PC:
```
eig.val <- snp.pca2$eig
eig.perc <- 100*snp.pca2$eig/sum(snp.pca2$eig)
eigen <- data.frame(eig.val,eig.perc)
write.csv(eigen,file="Qfly_adegenetPCA_eigen-summary.csv",row.names=TRUE,quote=FALSE)
```
Now, plot (taking the PC1 and PC2 percentages from the first two rows of the 'eig.perc' column in the Qfly_adegenetPCA_eigen-summary.csv file):
```
data2 <- read.delim("Qfly_adegenetPCA.txt")
mycol <- c("#f1c039","#f37d21", "#51692d", "#56ba32")
ggplot(data2, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(size = 2) + 
  theme_classic()+
  xlab("PC1 (1.99%)") +
  ylab("PC2 (1.34%)")
```
Now, let's recolour the points by status (first we need to redefine populations by creating a new genlight object based on pop.data2$V3):
```
gl.snp3 <- vcfR2genlight(snp_vcf2)
pop(gl.snp3) <- rep(pop.data2$V3)
```
Now repeat the steps from above to generate the PCA:
```
snp.pca3 <- glPca(gl.snp3, nf = 10)
snp.pca.scores3 <- as.data.frame(snp.pca3$scores)
snp.pca.scores3$pop <- pop(gl.snp3)
write.table(snp.pca.scores3, "Qfly_adegenetPCA_2.txt", sep = "\t")
eig.val3 <- snp.pca3$eig
eig.perc3 <- 100*snp.pca3$eig/sum(snp.pca3$eig)
eigen3 <- data.frame(eig.val3,eig.perc3)
write.csv(eigen3,file="Qfly_adegenetPCA_eigen-summary3.csv",row.names=TRUE,quote=FALSE)

data3 <- read.delim("Qfly_adegenetPCA_2.txt")
ggplot(data3, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(size = 2) + 
  theme_classic()+
  xlab("PC1 (1.99%)") +
  ylab("PC2 (1.34%)")
```
This PCA plot should show four main clusters, three of which are invasive (i.e., Alice Springs and the two Pacific Island populations); with the other corresponding to native and expanded populations in Australia.

Take a moment to think about whether this pattern meets your expectations. What does it say about how closely related the expanded and invasive populations are to the native ones? Does this make sense?

### Admixture
Let's look at admixture patterns using the LEA package in R. First, install the packages:
```
BiocManager::install("LEA", force = TRUE)
library(LEA)
library(remotes)
remotes::install_github('royfrancis/pophelper')
library(pophelper)
```
Note, you may first need to run:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
Now, create the input files for LEA's sNMF function:
```
vcf2geno(input.file = "Qfly.vcf", output.file = "Qfly.geno")
```
You should see three new files in your directory: Qfly.vcfsnp, Qfly.removed, Qfly.geno.
Do the clustering (see more here: https://www.rdocumentation.org/packages/LEA/versions/1.4.0/topics/snmf):
```
projectalpha <- NULL
projectalpha <- snmf("Qfly.geno", K = 1:10, repetitions = 50, entropy = TRUE, CPU = 8, project = "new")
```
Note, this is set to run on 8 threads (CPU = 8), which you might need to adjust for your system. 
Once this step is finished (it will take several minutes), you'll have a new folder in your directory ('Qfly.snmf'), which contains the output for all 10 K-values (with fifty runs per K-value).
The next step is to work out the cross-entropy criterion for all runs in the sNMF project:
```
plot(projectalpha, col = "maroon4", pch = 19, cex = 1.2)
```
This plot shows that 3 is the optimal number of ancestral populations, as that is the number of with the lowest cross-entropy value.
Now you can work out the best run (i.e., replicate) for each of the tested K-values:
```
best2 <- which.min(cross.entropy(projectalpha, K = 2))
best2
best3 <- which.min(cross.entropy(projectalpha, K = 3))
best3
```
You can do this all the way up to best10 by creating those variables (best4, best5, etc) following the same procedure used for best2 and best3.
Note that, due to variation in the algorithm, you may get different best runs from each other.
Now, let's plot the ancestry coefficients for our best run for K = 3:
```
my.col3 <- c("#51692d","#56ba32","#f1c039")
barchart(projectalpha, K = 3, run = best3,
         border = NA, space = 0,
         col = my.col3,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)
```
This plot shows us admixture profiles for all individuals. Let's colour it by status (i.e., native, expanded, invasive).
First, make the Q-matrix:
```
qmatrix <- as.data.frame(Q(projectalpha, K = 3, run = best3)) 
head(qmatrix)
```
Label column names of qmatrix:
```
ncol(qmatrix)
cluster_names <- c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) <- cluster_names
head(qmatrix)
```
Load in the population information:
```
popIds <- pop.data2
Ind <- popIds$V1
qmatrix$Ind <- Ind
Site <- popIds$V3
qmatrix$Site <- Site
```
Convert the dataframe to long format:
```
qlong <- melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)
```
Adjust the facet labels to show the three main groups:
```
levels(qlong$Site)
facet.labs <- c("Native", "Expanded", "Invasive")
levels(qlong$Site) <- facet.labs
levels(qlong$Site)
```
Define the colour palette:
```
pal <- colorRampPalette(c("#51692d","#56ba32","#f1c039")) 
cols <- pal(length(unique(qlong$variable)))
```
Order levels before plotting - you can change the order of these and it will change the order of how they are plotted:
```
qlong$Site <- ordered(qlong$Site, levels = c("Native", "Expanded", "Invasive")) 
qlong$Site
```
And plot:
```
ggplot(data=qlong, aes(x=Ind, y=value, fill = variable)) +
  geom_bar(stat="identity")+
  geom_col(color = "gray", size = 0.1) +
  scale_fill_manual(values = my.col3)+
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Admixture proportion") +
  facet_grid(~Site, scales="free", space="free")+
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(colour="black", size=9),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "none")
```
Examine the plot. What does it tell you about the evolutionary relationships and population structure among the three regions? Is it consistent with the other results so far?

Use the above code as a template and replot the results for a few different K-values, e.g., K=4 and K=6, to see how patterns might change.

### Population assignment
Finally, let's move into the adgenet package in R, using DAPC for population assignment. Note that DAPC can also be used to find clusters and generate PCAs, but we will use it just for population assignment here.
The first step is to create a training and 'supplementary' dataset from the genlight object. Here, we will randomly assign individuals to each category, but a very useful approach would be to train the algorithm on known samples and then use unknown invasive samples as the 'supplementary' individuals, to see if we can assign an origin for them with confidence (see this article for an example that applies population assignment using assignPOP for pest samples intercepted at the New Zealand border: https://pmc.ncbi.nlm.nih.gov/articles/PMC10099481/).
```
set.seed(2)
kept.id <- unlist(tapply(1:nInd(gl.snp2), pop(gl.snp2),
                         function(e) sample(e, 25,replace=T)))
x <- gl.snp2[kept.id]
x.sup <- gl.snp2[-kept.id]
nInd(x)
nInd(x.sup)
```
Now, run the DAPC analysis. Note that I did some initial steps to work out the best number of principal components (n.pca) and discriminant axes (n.da) to use here; see the manual for more information on how to run these steps (https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-dapc.pdf).
```
dapc1 <- dapc(x,n.pca=11,n.da=15)
pred.sup <- predict.dapc(dapc1, newdata=x.sup)
```
Plot the main data:
```
col <- funky(length(levels(pop(x))))
col.points <- transp(col[as.integer(pop(x))],.2)
scatter(dapc1, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, xlim=c(-10,10), legend=TRUE)
```
Plot the supplementary individuals on top:
```
par(xpd=TRUE)
points(dapc1$ind.coord[,1], dapc1$ind.coord[,2], pch=20,
       col=col.points, cex=3)
col.sup <- col[as.integer(pop(x.sup))]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=15,
       col=transp(col.sup,.7), cex=2)
```
See the percentage of individuals that were correctly assigned:
```
mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
```
The algorithm hasn't worked overly well. Only about 53% of individuals are correctly assigned.
Visualise the output:
```
table.value(table(pred.sup$assign, pop(x.sup)), col.lab=levels(pop(x.sup)))
```
Columns correspond to actual groups (origins) of supplementary individuals, while rows are inferred clusters.
Try re-running the above using the dataset where populations are defined by status. For example, here are the first steps:
```
set.seed(3)
kept.id2 <- unlist(tapply(1:nInd(gl.snp3), pop(gl.snp3),
                         function(e) sample(e, 25,replace=T)))
x2 <- gl.snp3[kept.id2]
x2.sup <- gl.snp3[-kept.id2]
nInd(x2)
nInd(x2.sup)
dapc2 <- dapc(x2,n.pca=11,n.da=15)
```
You should find an improvement in the percentage of corrrectly assigned individuals (to around 65%). You should also find a pattern in which the source for expanded populations is the native populations. The invaded populations show assignment to themselves, but we are looking at <25 supplementary individuals, so this is not surprising. As noted above, the best way to use this method in the context of inavsion genomics is to try to trace the source populations for intercepted border samples, or other invasive populations of unknown origin.

Congratulations - this is the official end of Day One! You should now have a good feel for how to generate FST, PCA, Admixture, and Population assignment plots for your own data, including the use of different population labelling techniques to best interpret your results.

### Bonus options for the extra keen:
1. Try to work out how to run the DAPC population assignment using only a few invasive individuals as 'supplementary individuals'. This will involve playing around with the x <- gl.snp2[kept.id2] part of the code.
2. Explore the DAPC manual and work out how to generate a PCA plot and a compoplot - how do they compare to the PCA/admixture methods used today?
3. Explore the R package assignPOP, which extends DAPC for population assignment using a machine-learning framework. For this, you would first have to use a program like VCFTOOLS to extract separate VCF files for each population of interest. You'd then need to use something like vcf2genepop.pl (from: https://github.com/z0on/2bRAD_denovo/blob/master/vcf2genepop.pl) to convert each VCF file into a genepop file. Finally, you'd need to run the assignment in R with assignPOP. 
There are many more steps, so this method will take some time to explore properly. See the tutorial here: https://alexkychen.github.io/assignPOP/

### Day One resources
https://www.nature.com/articles/s41437-023-00657-y

http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf

http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf

https://github.com/thibautjombart/adegenet/tree/master

https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-dapc.pdf
