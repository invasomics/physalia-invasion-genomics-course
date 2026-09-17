### Day Two
Today, we'll be working again with the Qfly dataset. We'll use two key methods to identify potential outliers that may be under selection in the different populations:

1. PCAdapt
2. BayPass

Note that we will have to run some of the steps outside of R, as BayPass is a standalone program.
   
Let's get into it!

### Load in the required libraries
First, set your working directory to something sensible, like 'DayTwo', then load in the required libraries.
Note, you may have to install them first with either install.packages(x) or BiocManager::install("x")!
```
library("pcadapt")
library("ggplot2")
library(qvalue)
library(dplyr)
library(UpSetR)
```
Read in the VCF file:
```
Qfly_vcf <- "Qfly.vcf"
```
Run PCadapt:
```
Qfly_pcadapt <- read.pcadapt(Qfly_vcf, type = "vcf")
```
Ignore the error about converting from vcf to pcadapt - it does still work currently; in future you might need to first convert your VCF file into BED file format (e.g., using PLINK software).
Now, run the pcadapt function and plot the resulting scree plot:
```
Qfly_pcadapt_kplot <- pcadapt(input = Qfly_pcadapt, K = 20)
plot(Qfly_pcadapt_kplot, option = "screeplot")
```
This plot is set to test 20 values of K and the resulting plot displays the percentage of variance explained by each principal component (PC) in decreasing order to help you choose the optimal number of components (K). You're looking for the point on the plot after which the curve starts to flatten out, so in this plot, a K-value of 3 is most appropriate. Thus, run a PCA using K=3:
```
Qfly_pcadapt_pca <- pcadapt(Qfly_pcadapt, K = 3)
```
Now check to see whether there is actually additional population structure on other PCs. First, read in the population data, then plot PCAs coloured by population:
```
pop.data2 <- read.table("Qfly_popmap.txt", header = F)
print(pop.data2)
plot(Qfly_pcadapt_kplot, option = "scores", i = 1, j = 2, pop = pop.data2$V2)
```
The PCA looks like we capture the population structure well with K=3 (see also results from Day One).
Now, let's plot a qq-plot and a manhattan plot:
```
plot(Qfly_pcadapt_pca, option = "qqplot")
plot(Qfly_pcadapt_pca, option = "manhattan")
```
The qqplot shows whether the dataset follows the expected uniform distribution under a neutral model. In this case, the plot confirms that most of the p-values follow the expected uniform distribution (i.e., the red line). However, the right hand side of the plot shows points that depart substantially from the diagonal red line. These represent outlier loci potentially under selection.

The manhattan plot shows the distribution of p-values for all SNPs that have a minor allele frequency > 0.05. The x-axis shows the SNPs, while the y-axis shows the negative base-10 logarithm of the p-values measuring the association between each SNP and the principal components.

Next, let's examine the p-value frequency data:
```
Qfly_pcadapt_pvalues <- as.data.frame(Qfly_pcadapt_pca$pvalues)
hist(Qfly_pcadapt_pca$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "grey")
```
This shows you the distribution of the p-values. Under the null hypothesis (neutral evolution), p-values should show a uniform distribution - in other words, the histogram bars should look relatively flat and equal across most of the chart (from 0.1 to 1.0). If there are loci under potential selection, you will see a spike in the histogram near 0.0. This indicates an excess of loci with very small p-values, which are your candidate outliers (markers under selection). 

Next, let's perform a Bonferroni correction for alpha=0.1:
```
Qfly_pcadapt_padj <- p.adjust(Qfly_pcadapt_pca$pvalues,method="bonferroni")
alpha <- 0.05
outliers <- which(Qfly_pcadapt_padj < alpha)
length(outliers)
```
This should identify 18 outliers at an alpha of 0.05 (and 20 if we set alpha to 0.10). View them with:
```
outliers
```
You can also print them to file:
```
write.table(outliers, file="Qfly_pcadapt_outliers.txt")
```
That's it for PCAdapt! We'll come back to the results later to see how they compare to BayPass.

### BayPass
Now, you need to move outside of R. The first step is to install BayPass. To do this, follow the instructions here for your system: https://forge.inrae.fr/mathieu.gautier/baypass_public.
Here's an example that worked on my linux machine, typed in a terminal window:
```
git clone https://forge.inrae.fr/mathieu.gautier/baypass_public.git
cd baypass_public/sources
make clean all FC=gfortran
```
Now, copy the resulting g_baypass file from the baypass_public/sources directory to your working directory. Note, your version of BayPass  may be named differently (e.g., ifx_baypass or i_baypass).

Now you're ready to run BayPass. 
Note that I've provided the Qfly.bp (which I created using: python ./reshaper_baypass.py Qfly.vcf Qfly_popmap.txt Qfly.bp; the python script is available here: https://gitlab.com/YDorant/Toolbox/-/tree/master?ref_type=heads).
To run BayPass you also need to have specific pairwise constrast files set up. This first example shows how to run BayPass using all native populations vs the invasive Alice population. You set this contrast up in BayPass by using a '1' to indicate the native populations, a '0' to indicate populations that should be ignored in the current analysis, and a '-1' to indicate the population for contrast (in this case, Alice). The ecotypes files you will need are all provided in the DayTwo folder.
In a terminal window, type:
```
./g_baypass -gfile Qfly.bp -contrastfile Qfly_NativevsAlice.ecotype -efile Qfly_NativevsAlice.ecotype -outprefix Qfly_NativevsAlice -nthreads 6 -seed 2026
```
Note that you might need to change the number of threads to suit your computer (i.e., adjust -nthreads 6). 
This analysis will take anywhere from ~2 mins to 30 mins or more, depending on the specs of your computer.
This analysis will produce a series of 9 output files with the prefix 'Qfly_NativevsAlice'. Of these, we are most interested in the 'Qfly_NativevsAlice_summary_contrast.out' file, which we will now read into R.
If for some reason you cannot get BayPass to install and work, I have provided the three key BayPass output files in the DayTwo/RequiredFiles/BackupFiles directory. But, please do try to get BayPass to work on your own!

#### Return to R
Read in the contrast file and rename column 4 to be a bit more useful:
```
Alice.C2 <- read.table("Qfly_NativevsAlice_summary_contrast.out",h=T) 
colnames(Alice.C2)[4]  <- "Log10pvals"
```
Now, import the provided (in the DayTwo/RequiredFiles/ directory) scaffold list, which tells us where the loci are located on the various scaffolds. This file was created by running some commands in bash on the Qfly.vcf file: cat Qfly.vcf | grep -v "#" | cut -f3 > Qfly_scaffold_list.txt
Import this file into R, and append the column to the Alice.C2 dataframe:
```
scaffolds <- read.table("Qfly_scaffold_list.txt")
Alice.C2 <- as.data.frame(cbind(Alice.C2, scaffolds))
colnames(Alice.C2)[5] <- "Scaffold"
```
Now, let's check the behaviour of the p-values associated with the C2 statistic:
```
hist(10**(-1*Alice.C2$Log10pvals),freq=F,breaks=40)
abline(h=1)
```
This plot is similar to the histogram we created above for PCAdapt, which we later followed with Bonferroni correction. Here, we will convert the p-values from log10 to normal, and then convert to q-values:
```
pvalues <- as.data.frame(10^-Alice.C2[,4]) 
Alice.C2 <- cbind(Alice.C2,pvalues) 
colnames(Alice.C2)[6]  <- "pvals"
pvalues <- as.vector(Alice.C2$pvals)
qval <- qvalue(p = pvalues)
plot(qval$qvalues)
abline(h=0.05)
```
Q-values are adjusted p-values used to control the False Discovery Rate (FDR) when identifying outlier loci under selection. Q-values < 0.05 are highly significant. Let's extract the SNPs that have q-values < 0.05 and write them to a file:
```
qvalues <- as.data.frame(qval$qvalues)
Alice.C2 <- cbind(Alice.C2,qvalues)
colnames(Alice.C2)[7]  <- "Qvals"

selected_SNPs <- Alice.C2[Alice.C2$Qvals < 0.05, ]
write.table(selected_SNPs,"BPoutliers_NativevsAlice_FDR5.txt", sep = "\t")
```
There should be ~57 outliers. Though different BayPass runs may vary in the absolute number of potential outliers detected because of stochasticity in the algorithm or use of different software versions/computer platforms, you should see a similar trend in your plots to those in the manuscript. Let's create the outlier plot:
```
ggplot(Alice.C2, aes(x=MRK, y = Log10pvals)) + 
  geom_point(show.legend = FALSE, alpha = 1, size = 2) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  theme_classic()
```
We can also replot these by chromosome, using the provided file, where scaffold names are updated to chromosome numbers (in the file, the number 6 is used for unplaced scaffolds):
```
scaffolds2 <- read.table("Qfly_chromosomes-chrNumbers.txt", header = TRUE) 
Alice.C2_2 <- cbind(Alice.C2, scaffolds2)
ggplot(Alice.C2_2, aes(x=MRK, y = Log10pvals, color = CHR)) + 
  geom_point(show.legend = FALSE, alpha = 1, size = 2) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  theme_classic()
```

#### Next steps
Re-run the above steps for the other two BayPass comparisons. In BayPass (outside of R), run:
```
./g_baypass -gfile Qfly.bp -contrastfile Qfly_NativevsExpanded.ecotype -efile Qfly_NativevsExpanded.ecotype -outprefix Qfly_NativevsExpanded -nthreads 6 -seed 2026
./g_baypass -gfile Qfly.bp -contrastfile Qfly_NativevsIslands.ecotype -efile Qfly_NativevsIslands.ecotype -outprefix Qfly_NativevsIslands -nthreads 6 -seed 2026
```
Then follow the rest of the code above to generate outlier lists for these two additional comparisons, creating these two BayPass output files:
1. Qfly_NativevsExpanded_summary_contrast.out
2. Qfly_NativevsIslands_summary_contrast.out
   
Then, back in R, repeat the rest of the code to generate the outlier files:
1. BPoutliers_NativevsExpanded_FDR5.txt
2. BPoutliers_NativevsIslands_FDR5.txt

#### Common outliers
Now, make a Venn diagram to check for common SNPs between the three BayPass comparisons:
```
Native <- read.table("BPoutliers_NativevsAlice_FDR5.txt", header=T)
Expanded <- read.table("BPoutliers_NativevsExpanded_FDR5.txt", header=T)
Islands <- read.table("BPoutliers_NativevsIslands_FDR5.txt", header=T)

library(ggvenn)
venn1 <- list(Native=Native$MRK, Expanded=Expanded$MRK, Islands=Islands$MRK)
ggvenn(venn1, fill_color = c("#F8D210", "#593F1F", "#DEE3CA"), fill_alpha = 0.7, stroke_size = 0.2, set_name_size = 4, stroke_color = "black", show_percentage = FALSE)
```

Let's also check to see how the BayPass result compares with PCAdapt. We already generated the table of PCAdapt outliers above, so:
```
PCadapt <- read.table("Qfly_pcadapt_outliers.txt")
venn2 <- list(Native=Native$MRK, Expanded=Expanded$MRK, Islands=Islands$MRK, PCAdapt=PCadapt$x)
ggvenn(venn2, fill_color = c("#F8D210", "#593F1F", "#DEE3CA", "grey"), fill_alpha = 0.7, stroke_size = 0.2, set_name_size = 4, stroke_color = "black", show_percentage = FALSE)
```
Your results should look similar to those provided in the DayTwo folder (BayPassVenn.jpeg and FullVenn.jpeg), with some potential absolute differences driven by computer platform and algorithm run variation. You'll see that the results are a little different to those in the published article (which only compares the BayPass results among the different introduced clusters and also ran a more extensive BayPass analysis overall).
Another nice way to plot these results is here:
```
upset(
  data = fromList(venn2), 
  order.by = "freq", 
  empty.intersections = "on", 
  point.size = 3.5, 
  line.size = 2, 
  mainbar.y.label = "Outlier Count", 
  sets.x.label = "Total Outliers", 
  text.scale = c(1.3, 1.3, 1, 1, 2, 1.3), 
  number.angles = 30, 
  nintersects = 11
) 
```

A final note on BayPass is that we've only run it once per contrast today. For robust use, it needs to be run 3-5 times per tested contrast, each time with a different seed, then the values in the output omega file should be checked to make sure the results are consistent.

Congratulations - this is the official end of Day Two! You should now have a good feel for how to generate information on candidate adaptive outliers for your own data, using two different methods (PCAdapt and BayPass).

### A note on outlier annotation
Functional annotation of candidate outlier genes is a big task and we don't have time to complete those steps today. For more information on how to complete those steps if you have genome annotation files available for your own data, see the tutorial here: https://github.com/Elahep/B.tryoni_PopGenomics/tree/main/4-GOterms

### Bonus steps
1. Try running PCAdapt with K=4 instead of K=3. Do the results change very much?
2. If you wish to compare your results more directly with the Qfly publication (https://www.nature.com/articles/s41437-023-00657-y#MOESM1), try also generating BayPass results for a 'Native vs All invasive' population contrast. To achieve this, you will have to create an ecotype file, placing 1, 0, and -1 so that you have all Native populations coded as '1', all Expanded populations coded as '0', and the three invasive populations coded as '-1'. See if you can replicate the Venn diagram in Figure 3b of the manuscript (though your exact numbers may vary a litte, the relative patterns should be consistent). 

### Resources for today:
https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

https://forge.inrae.fr/mathieu.gautier/baypass_public
