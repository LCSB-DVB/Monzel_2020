###############################################################################
## (c) T.KAOMA, BioMod, Luxembourg Institute of Health, 2019-08-27
## TOX-Project: Analyze of data. Analyze included: 
##     > log transformation
##     > Normalization: z-score per feature within each expr/CL/section
##     > unsupervised analysis & heatmap
##     > Classification with Random Forest
##     > Export of figures and results
###############################################################################


rm(list=ls())
graphics.off()

suppressMessages({
    library(tidyr)
    library(tibble)
    library(dplyr)
    library(tidyverse)
    library(RColorBrewer)
    library(pheatmap)
    library(randomForest)
    library(parallel)
    library(gridExtra)
    library(grid)
    library(ComplexHeatmap)
})


file_idx="2.2-ECS-"
wd = "~/projects/ML-Image-LCSB/2019-08-26/"
path2data = "data/20190408_Data_ALL_Cleaned_final.csv"
script = wd

stopifnot(dir.exists(wd))
stopifnot(file.exists(path2data))
stopifnot(file.exists(file.path(script,"1-pvca.R")))

source(file.path(script,"1-pvca.R"))

setwd(wd)
empty.plot <- function() plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes = F,xlab="",ylab="")

#####
Data = read.csv(path2data,as.is = T,check.names = F)
head(Data)
metadata = Data[,1:which(colnames(Data) == "Section")]
rownames(metadata) = metadata$ID
exprs = Data[,(which(colnames(Data) == "Section")+1):ncol(Data)]
rownames(exprs) = metadata$ID
exprs = exprs[,order(colnames(exprs))]

#### Only CTRL and 60HDA
metadata = metadata[metadata$Condition %in% c("CTRL","6OHDA"),]
exprs    = log(exprs[rownames(metadata),])
range(exprs)
dim(exprs)

exprs_norm_split = split(1:nrow(metadata),with(metadata,paste(Experiment,Cellline,Section)))
exprs_norm_split = exprs_norm_split[sapply(exprs_norm_split,length) > 1]
exprs_norm = do.call(rbind,lapply(
    exprs_norm_split, function (x){
        exprs_ = exprs[x,]
       
         ### overall mean
        oa_mean = apply(exprs_,2,mean)
        oa_sd   = apply(exprs_,2,sd)
        
        ### zscore by feature
        zscore = sweep(exprs_,2,oa_mean,"-")
        zscore = sweep(zscore,2,oa_sd,"/")
        
        ### merge meta and zscore
        zscore
    }
))
exprs = exprs_norm
metadata = do.call(rbind,lapply(
    exprs_norm_split, function (x){
        metadata[x,]
    }
))
all(rownames(exprs) == rownames(metadata))
rm(exprs_norm)

###########
spl.cols = c(
    CTRL     = "lightblue",
    "6OHDA"  = "red"
)[metadata$Condition]

feature.names = c(
    MAP2ByNuc='MAP2+/Nuclei',
    MAP2ByNucAlive ='MAP2+/Nuclei alive',
    MAP2MaskSum    = 'MAP2 mask',
    NucAliveByNuc  = 'Nuclei alive/Nuclei',
    NucDeadByNuc   = 'Nuclei dead/Nuclei',
    NucMaskAlive   = 'Nuclei alive mask',
    NucMaskHigh    = 'Nuclei dead mask',
    NucMaskSum     = 'Nuclei mask',
    SkelTH         = 'TH skeleton',
    THByMap2       = 'TH+/MAP2+',
    THByNuc        = 'TH+/Nuclei',
    THByNucAlive   = 'TH+/Nuclei alive',
    THByTuj1       = 'TH+/TUJ1+',
    THFragmentation='TH fragmentation',
    THMaskSum      = 'TH mask',
    THPercent      = 'Percentage TH',
    Tuj1ByNuc      = 'TUJ1+/Nuclei',
    Tuj1ByNucAlive = 'TUJ1+/Nuclei alive',
    Tuj1MaskSum    = 'TUJ1 mask'
)

exprs = exprs[,names(feature.names)]

fact.cols = list(
    Section   = c("Border" = "gold1", "Center" = "gold4"),
    Condition = c(CTRL="lightblue","6OHDA"="red"),
    Cellline  = c("34769"="darkgrey",K7="grey80",T129="grey50")
)

########### Quality control and data exploration
explore.file = file.path("Results",paste0(file_idx,"QC.pdf"))

pdf(explore.file,width = 11.4,height = 8.4)
par(mar=c(10,6,6,6))
boxplot(
    t(exprs),
    col=spl.cols,
    las=2,
    main="Boxplot of observations - log(Norm - ESC)",
    cex.axis=.5
)

dens = apply(exprs,1,density)
dx = sapply(dens,function(x)x$x)
dy = sapply(dens,function(x)x$y)
matplot(
    dx,dy,
    col=spl.cols,
    lty=1,
    type='l',
    lwd=1,
    main="Density plot of observations - log(Raw data)"
)
legend(
    "topright",
    legend = names(fact.cols$Condition),
    fill   = fact.cols$Condition,
    bty="n"
)

########## All features - PVCA
batch.factors = c("Experiment","Cellline","Replicate","Condition","Section")
cmb = NULL
if(length(batch.factors) > 1) {
    cmb = combn(batch.factors,2)
    cmb = paste(cmb[1,],cmb[2,],sep=":")
    cmb = cmb[!grepl("Replicate",cmb)]
}
pct_threshold <- 0.9
model.func = paste("(1|", c(batch.factors,cmb), ")", sep = "")
pvcaObj <- pvca(
    as.matrix(t(exprs)),
    metadata,
    batch.factors,
    pct_threshold,
    model.func
)
fac.order = c(batch.factors,cmb,"resid")

par(mar=c(10,6,6,6))
bp <- barplot(
    pvcaObj$dat[match(fac.order,pvcaObj$label)], xlab = "",
    ylab = "Weighted average proportion variance",
    ylim= c(0,1.1),col = c("blue"), las=2,
    main="PVCA estimation bar chart"
)
text(bp,-0.02, srt = 45, adj= 1, xpd = TRUE, labels = fac.order, cex=0.7)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat[match(fac.order,pvcaObj$label)],labels = new_values[match(fac.order,pvcaObj$label)], pos=3, cex = .8)

pct_threshold <- 0.9
model.func = paste("(1|", c(batch.factors), ")", sep = "")
pvcaObj <- myPVCA(t(exprs)[,], metadata[,], batch.factors, pct_threshold, model.func)
fac.order = c(batch.factors,"resid")

par(mar=c(10,6,6,6))
bp <- barplot(
    pvcaObj$dat[match(fac.order,pvcaObj$label)], xlab = "",
    ylab = "Weighted average proportion variance",
    ylim= c(0,1.1),col = c("blue"), las=2,
    main="PVCA estimation bar chart"
)
text(bp,-0.02, srt = 45, adj= 1, xpd = TRUE, labels = fac.order, cex=0.7)
values = pvcaObj$dat
new_values = round(values , 3)
text(
    bp,
    pvcaObj$dat[match(fac.order,pvcaObj$label)],
    labels = new_values[match(fac.order,pvcaObj$label)],
    pos=3,
    cex = .8
)

######
spl.cor = cor(t(exprs))
fHclust = hclust(as.dist(1-spl.cor))
sHclust = hclust(as.dist(1-spl.cor))

Heatmap(
    spl.cor,
    cluster_rows = fHclust,
    cluster_columns  = sHclust,
    top_annotation = HeatmapAnnotation(
        metadata[,c("Section","Condition","Cellline","Experiment")],
        col = fact.cols,
        show_annotation_name = T
    ),
    show_column_names = F,
    show_row_names = F
)


spl.cor = cor(t(scale(exprs)))
fHclust = hclust(as.dist(1-spl.cor))
sHclust = hclust(as.dist(1-spl.cor))

Heatmap(
    spl.cor,
    cluster_rows = fHclust,
    cluster_columns  = sHclust,
    top_annotation = HeatmapAnnotation(
        metadata[,c("Section","Condition","Cellline","Experiment")],
        col = fact.cols,
        show_annotation_name = T
    ),
    show_column_names = F,
    show_row_names = F
)

dev.off()
Biobase::openPDF(explore.file)

########### B. Random Forest/ML
factor.names = c("Experiment","Cellline","Replicate","Condition","Section")

#### 0. Requirement
set.seed(123)
nrows = nrow(exprs)
nfolds = 5 # 5-fold CV
cv_times = 10 # 10 x 5-fold CV

table(metadata$Condition)
# 6OHDA  CTRL 
# 52    46
folds_ll = lapply(1:cv_times,function(zz) lapply(
    split(1:nrows,metadata$Condition),
    function(x){
        split(x,sample(cut_number(x,n=nfolds,label=F)))
    }
))

lapply(folds_ll,lapply,sapply,length)
foldss = lapply(folds_ll,function(folds_l) mapply(
    function(a,b){
        sort(unlist(
            c(folds_l[[1]][a],folds_l[[2]][b])
        ))
    },
    sample(1:5),
    sample(1:5),
    SIMPLIFY = F
))

pdf.output = file.path("Results",paste0(file_idx,"ML-RF.pdf"))
pdf(pdf.output,width = 13,height = 8.4)

##### 1. Condition
empty.plot()
text(.5,.5,"******* Condition *******",cex=2)

RF = list()
### 1. Level per fold"
for(idx in 1:length(foldss)){
    folds = foldss[[idx]]
    cat("Run ",idx,"/10\n")
    
    ## Remove HDA and NAC as factors
    dF = data.frame(Condition = metadata[,"Condition"],exprs)
    dF[!sapply(dF, is.numeric)] = lapply(dF[!sapply(dF, is.numeric)],as.factor)
    
    ## Prediction with RandomForest
    rF.list = mclapply(1:length(folds),function(i){
        tx   = sort(unlist(folds[-i]))
        test = folds[[i]]
        rF   = randomForest(Condition ~., data=dF[tx,])
        pred = predict(rF,dF[test,])
        list(
            rF = rF,
            accuracy   = mean(dF$Condition[test] == pred),
            importance = rF$importance[order(rF$importance, decreasing = T),],
            confusion  = table(Obs = dF$Condition[test], pred = pred)
        )
    },
    mc.cores= detectCores() - 1
    )
    RF[[idx]] = rF.list
    rm(folds)
}

### Summary
mat = round(do.call(rbind, lapply(RF, sapply,function(x) x$accuracy)),2)
colnames(mat) = paste0('Fold ',1:nfolds)
rownames(mat) = paste0('Run ',1:cv_times)

par(mfrow=c(1,1))
empty.plot()
text(
    .5,.97,
    paste0(
        "Condition"
    ),cex=2
)

text(
    .5,.90,
    paste0(
        "Random Forest (10x 5-fold CV)"
    ),cex=1.2
)

text(
    .5,.86,
    paste0(
        "\nMean of accuracy = ",round(mean(rowMeans(mat)),2),
        " ± ",
        round(sd(rowMeans(mat)),2),"\n"
    ),cex=1.2
)
grid.table(mat)

layout(mat = matrix(c(1,2)),heights = c(.2,.8))
par(mar=c(4,6,2,6),xpd=T)
# empty.plot()
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes = F,xlab="",ylab="")
text(
    .5,.005,
    paste0(
        "Distribution of RF accuracies\nover 10x 5-fold CV\n",
        "Grand mean = ",round(mean(rowMeans(mat)),2),
        " ± ",
        round(sd(rowMeans(mat)),2),"\n"
    ),cex=1.2
)
boxplot(lapply(RF, function(x) sapply(x,function(y) y$accuracy)))

df = data.frame(
    Run = rownames(mat),
    Accuracy = rowMeans(mat),
    sd = apply(mat,1,sd)
)
par(mfrow=c(1,1),mar=c(6,6,6,4))
x = barplot(
    df$Accuracy,
    ylim=c(0,max(df$Accuracy + df$sd)+.05),
    ylab = "Average accuracy",
    main=paste0(
        "\n\nAverage accuracy\n(10x 5-fold CV)\n",
        "Grand mean = ",round(mean(rowMeans(mat)),2),
        " ± ",
        round(sd(rowMeans(mat)),2),"\n"
    ),main.cex=1.2
)
for(i in 1:nrow(df)) lines(c(x[i],x[i]),with(df,c(Accuracy[i]-sd[i],Accuracy[i]+sd[i])))
for(i in 1:nrow(df)) lines(c(x[i]-.1,x[i] + .1),rep(with(df,c(Accuracy[i]-sd[i])),2))
for(i in 1:nrow(df)) lines(c(x[i]-.1,x[i] + .1),rep(with(df,c(Accuracy[i]+sd[i])),2))
axis(1,x,df$Run,las=2)
rm(x,df)

vars = names(RF[[1]][[1]]$importance)
idx = 0
imports = list()
for(rF.list in RF)
    imports[[idx <<- idx + 1]] = sapply(rF.list, function(x) x$importance[vars])
imports = do.call(cbind, lapply(imports, function(x) x))

rownames(imports)  = feature.names[rownames(imports)]


par(mfrow=c(1,1),mar=c(4,10,6,4))
barplot(
    sort(rowMeans(imports),decreasing = F),
    main="\nRF - Condition - Feature importances",
    horiz = T,las=2,xlab = "Feature importances\n"
)
model1 = list(
    RF = RF,
    imports = imports
)
rm(mat,RF,imports,vars,mat,idx)

##### ##### 2. Condition - with ECS
empty.plot()
text(.5,.5,"******* Condition (features + ECS) *******",cex=2)

RF = list()
### 1. Level per fold"
for(idx in 1:length(foldss)){
    folds = foldss[[idx]]
    cat("Run ",idx,"/10\n")
    
    ## Remove HDA and NAC as factors
    dF = data.frame(metadata[,c("Experiment","Cellline","Section","Condition")],exprs)
    dF[!sapply(dF, is.numeric)] = lapply(dF[!sapply(dF, is.numeric)],as.factor)
    
    ## Prediction with RandomForest
    rF.list = mclapply(1:length(folds),function(i){
        tx   = sort(unlist(folds[-i]))
        test = folds[[i]]
        rF   = randomForest(Condition ~., data=dF[tx,])
        pred = predict(rF,dF[test,])
        list(
            rF = rF,
            accuracy   = mean(dF$Condition[test] == pred),
            importance = rF$importance[order(rF$importance, decreasing = T),],
            confusion = table(Obs = dF$Condition[test], pred = pred)
        )
    },
    mc.cores= detectCores() - 1
    )
    RF[[idx]] = rF.list
    rm(folds)
}
rm(rF.list)

### Summary
mat = round(do.call(rbind, lapply(RF, sapply,function(x) x$accuracy)),2)
colnames(mat) = paste0('Fold ',1:nfolds)
rownames(mat) = paste0('Run ',1:cv_times)

par(mfrow=c(1,1))
empty.plot()
text(
    .5,.97,
    paste0(
        "Condition"
    ),cex=2
)

text(
    .5,.90,
    paste0(
        "Random Forest (10x 5-fold CV)"
    ),cex=1.2
)

text(
    .5,.86,
    paste0(
        "\nMean of accuracy = ",round(mean(rowMeans(mat)),2),
        " ± ",
        round(sd(rowMeans(mat)),2),"\n"
    ),cex=1.2
)
grid.table(mat)

layout(mat = matrix(c(1,2)),heights = c(.2,.8))
par(mar=c(4,6,2,6),xpd=T)
# empty.plot()
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes = F,xlab="",ylab="")
text(
    .5,.005,
    paste0(
        "Distribution of RF accuracies\nover 10x 5-fold CV\n",
        "Grand mean = ",round(mean(rowMeans(mat)),2),
        " ± ",
        round(sd(rowMeans(mat)),2),"\n"
    ),cex=1.2
)
boxplot(lapply(RF, function(x) sapply(x,function(y) y$accuracy)))

df = data.frame(
    Run = rownames(mat),
    Accuracy = rowMeans(mat),
    sd = apply(mat,1,sd)
)
par(mfrow=c(1,1),mar=c(6,6,6,4))
x = barplot(
    df$Accuracy,
    ylim=c(0,max(df$Accuracy + df$sd)+.05),
    ylab = "Average accuracy",
    main=paste0(
        "\n\nAverage accuracy\n(10x 5-fold CV)\n",
        "Grand mean = ",round(mean(rowMeans(mat)),2),
        " ± ",
        round(sd(rowMeans(mat)),2),"\n"
    ),main.cex=1.2
)
for(i in 1:nrow(df)) lines(c(x[i],x[i]),with(df,c(Accuracy[i]-sd[i],Accuracy[i]+sd[i])))
for(i in 1:nrow(df)) lines(c(x[i]-.1,x[i] + .1),rep(with(df,c(Accuracy[i]-sd[i])),2))
for(i in 1:nrow(df)) lines(c(x[i]-.1,x[i] + .1),rep(with(df,c(Accuracy[i]+sd[i])),2))
axis(1,x,df$Run,las=2)
rm(x,df)

vars = names(RF[[1]][[1]]$importance)
idx = 0
imports = list()
for(rF.list in RF)
    imports[[idx <<- idx + 1]] = sapply(rF.list, function(x) x$importance[vars])
imports = do.call(cbind, lapply(imports, function(x) x))

rownames(imports)[rownames(imports) %in% names(feature.names)]  = 
    feature.names[rownames(imports)[rownames(imports) %in% names(feature.names)]]

par(mfrow=c(1,1),mar=c(4,10,6,4))
barplot(
    sort(rowMeans(imports),decreasing = F),
    main="\nRF - Condition - Feature importances",
    horiz = T,las=2,xlab = "Feature importances\n"
)
model2 = list(
    RF = RF,
    imports = imports
)
rm(mat,RF,imports,vars,mat,idx)

dev.off()
Biobase::openPDF(pdf.output)

save(
    nfolds,cv_times,folds_ll,foldss,
    model1, model2,
    file=file.path("Results",paste0(file_idx,"ML-RF.RData"))
)