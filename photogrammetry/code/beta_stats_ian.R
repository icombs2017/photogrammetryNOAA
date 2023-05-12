library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')
library("usedist")
library("ggplot2")
library("microbiome")
library("multcompView")

rm(list=ls())
setwd("~/NFWF/090221EMillcus515F-295786525") 
mapfile = "~/NFWF/NFWF_metadata.txt"


pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="bonferroni",perm=999){
  
  library(vegan)
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    resp <- as.matrix(x)[sub_inds,sub_inds]
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

load(file = "ps_rename_rare.RData") #renamed, rarefied data
ps <- ps_rarefied
rm(ps_rarefied)
sample_sums(ps) #5000

mapfile = "~/NFWF/NFWF_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

ps <- subset_samples(ps, Coral_Species == "ACER")
nsamples(ps)

### Calculate distance matrices
# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
ps_bc <- phyloseq::distance(ps, method = "bray")
ps_bj <- distance(ps, method = "jaccard", binary =TRUE)

#PCoA plot
ps_bc <- phyloseq::distance(ps, method = "bray")


logt  = transform_sample_counts(ps, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues

plot_ordination(logt, out.pcoa.logt, type = "samples", 
                color = "Nutrient_no_level") + labs(col = "Nutrient") +
  coord_fixed(sqrt(evals[2] / evals[1]))

# Prepare data for plotting
meta <- as.data.frame(sample_data(ps_not_t0))
table <- as.data.frame(otu_table(ps_not_t0))
head(meta)

#these are for ellipses
ctrl <- rownames(meta[which(meta[,14] =="Ctrl"),])
T0 <- rownames(meta[which(meta[,14] =="T0"),])
nut <- rownames(meta[which(meta[,14] !="Ctrl" & meta[,14] !="T0"),])
ammonium <- rownames(meta[which(meta[,14] =="Ammonium"),])
nitrate <- rownames(meta[which(meta[,14] =="Nitrate"),])
phosphate <- rownames(meta[which(meta[,14] =="Phosphate"),])
combined <- rownames(meta[which(meta[,14] =="Combined"),])
three <- rownames(meta[which(meta[,6] =="3"),])
six <- rownames(meta[which(meta[,6] =="6"),])
Nutrient <- meta$Nutrient_no_level
Weeks <- meta$Exposure_weeks


# params for plotting
dims <- c(1,2)
ellp.kind <- "ehull"

# ordinate with Bray Curtis
object <- metaMDSiter(ps_bc, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE) #wow finally reached solution at run 977 stress 0.1109059 
save(object, file = "ord_bc_50.RData")

mds.fig <- ordiplot(object, display = "sites", type = "none", choices = dims)
plot(1, type="n", main= "NMDS All Timepoints", xlab="MDS1", ylab="MDS2", xlim=c(-4.5,1.4), ylim=c(-3,3))
points(mds.fig, "sites", pch = 18, col = "#1E88E5", select = T0)
points(mds.fig, "sites", pch = 17, col = "#D81B60", select = nut)
points(mds.fig, "sites", pch = 19, col = "#009E73", select = ctrl)
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#009E73", lwd = 2, show.groups = "Ctrl")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#D81B60", lwd = 2, show.groups = "Combined")
legend("bottomleft", legend = c("T0", "Control", "Nutrient"), pch = c(18, 19, 17), col = c("#1E88E5","#009E73","#D81B60"))
text(0.7,2.8, labels = c("stress=0.111"))

mds.fig <- ordiplot(object, display = "sites", type = "points", choices = dims)
plot(1, type="n", main= "NMDS 3 and 6 Weeks", xlab="MDS1", ylab="MDS2", xlim=c(-3.75,3), ylim=c(-3,2.2))
points(mds.fig, "sites", pch = 20, col = "#BE0032", select = ammonium)
points(mds.fig, "sites", pch = 19, col = "#F3C300", select = nitrate)
points(mds.fig, "sites", pch = 18, col = "#0067A5", select = phosphate)
points(mds.fig, "sites", pch = 15, col = "#009E73", select = combined) #square
points(mds.fig, "sites", pch = 17, col = "#848482", select = ctrl) #triangle
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#848482", lwd = 2, show.groups = "Ctrl")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#0067A5", lwd = 2, show.groups = "Phosphate")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#BE0032", lwd = 2, show.groups = "Ammonium")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#F3C300", lwd = 2, show.groups = "Nitrate")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#009E73", lwd = 2, show.groups = "Combined")
legend("bottomleft", legend = c("Control", "Ammonium", "Nitrate", "Phosphate", "Combined"), pch = c(17, 20, 19, 18, 15), col = c("#848482","#BE0032","#F3C300", "#0067A5", "#009E73"))
text(2.3,2, labels = c("stress=0.08"))

mds.fig <- ordiplot(object, display = "sites", type = "points", choices = dims)
plot(1, type="n", main= "NMDS 3 and 6 Weeks", xlab="MDS1", ylab="MDS2", xlim=c(-3.75,3), ylim=c(-3,2.2))
points(mds.fig, "sites", pch = 20, col = "#BE0032", select = three)
points(mds.fig, "sites", pch = 15, col = "#009E73", select = six)
ordiellipse(object, Weeks, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#BE0032", lwd = 2, show.groups = "3")
ordiellipse(object, Weeks, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#009E73", lwd = 2, show.groups = "6")
legend("bottomleft", legend = c("Three Weeks", "Six Weeks"), pch = c(20, 15), col = c("#BE0032", "#009E73"))
text(2.3,2, labels = c("stress=0.08"))
