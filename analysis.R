install.packages('meta')
library(meta)
library(ggplot2)
library(gridExtra)
library(grid)

########################################################
###             1. Read in datasets                  ###
########################################################

meta = read.csv('C:/Users/Ben/Desktop/Work/EBV review/EBV meta-analysis1910_seropositivity.csv')
meta_seropositives = read.csv('C:/Users/Ben/Desktop/Work/EBV review/EBV meta-analysis1910_seropositivity_pooled.csv')
meta_IM = read.csv('C:/Users/Ben/Desktop/Work/EBV review/EBV meta-analysis1910_IM.csv')
meta_IM = meta_IM[c(1:19),]
meta_PCR = read.csv('C:/Users/Ben/Desktop/Work/EBV review/EBV meta-analysis1910_PCR.csv')



########################################################
###             2. Meta-analyses                     ###
########################################################


#Meta for adults, serum, VCA

meta1 = metabin(data=meta,
                event.e=meta$nMS.,
                n.e=meta$nMS,
                event.c=meta$nCONTROL.,
                n.c=meta$nCONTROL,
                studlab = paste(meta$Author,meta$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5, 
                subset=(meta$Adults.paeds=='Adult' & meta$serum.CSF=='S' & meta$Ab == 'VCA')
                )
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/vcapositivity.bmp',res=300,width=8,height=8,units='in')
forest1 = forest(meta1, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
       )
dev.off()


funnel(meta1, 
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL,
       main = 'Funnel plot of VCA seropositivity'
       )


metabias(meta1)

#Meta for adults, serum, EBNA

meta2 = metabin(data=meta,
                event.e=meta$nMS.,
                n.e=meta$nMS,
                event.c=meta$nCONTROL.,
                n.c=meta$nCONTROL,
                studlab = paste(meta$Author,meta$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5, 
                subset=(meta$Adults.paeds=='Adult' & meta$serum.CSF=='S' & meta$Ab == 'EBNA1')
)

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/ebnapositivity.bmp',res=300,width=8,height=8,units='in')

forest(meta2, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)
forest2
dev.off()

funnel(meta2,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL,
       main = 'Funnel plot of EBNA seropositivity')

metabias(meta2)

#Meta for seropositive adults - nb where possible a metric like VCA+EBNA positivity has been used. If that was not available from the study, EBNA was used. 
meta4 = metabin(data=meta_seropositives,
                event.e=meta_seropositives$nMS.,
                n.e=meta_seropositives$nMS,
                event.c=meta_seropositives$nCONTROL.,
                n.c=meta_seropositives$nCONTROL,
                studlab = paste(meta_seropositives$Author,meta_seropositives$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5, 
                subset=(meta_seropositives$Adults.paeds=='Adult' & meta_seropositives$serum.CSF=='S')
)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/adult_seropositivity.bmp',res=300,width=8,height=8,units='in')

forest(meta4, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)
dev.off()


funnel(meta4,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL,
       main = 'Funnel plot of EBV seropositivity (adults)')

metabias(meta4)



#Meta for seropositive children - nb where possible a metric like VCA+EBNA positivity has been used. If that was not available from the study, EBNA was used. 

meta5 = metabin(data=meta_seropositives,
                event.e=meta_seropositives$nMS.,
                n.e=meta_seropositives$nMS,
                event.c=meta_seropositives$nCONTROL.,
                n.c=meta_seropositives$nCONTROL,
                studlab = paste(meta_seropositives$Author,meta_seropositives$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5, 
                subset=(meta_seropositives$Adults.paeds=='Paeds' & meta_seropositives$serum.CSF=='S')
)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/seropositivitykids.bmp',res=300,width=8,height=8,units='in')

forest(meta5, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)
dev.off()
grid::grid.text("Meta-analysis of EBV seropositivity (children)", 0.5, 0.95,gp=grid::gpar(fontsize=14))

funnel(meta5,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL)

metabias(meta5, k.min = 5)


#Meta for IM 
meta6 = metabin(data=meta_IM,
                event.e=meta_IM$nMS.,
                n.e=meta_IM$nMS,
                event.c=meta_IM$nControl.,
                n.c=meta_IM$nControl,
                studlab = paste(meta_IM$Author,meta_IM$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5)

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/IM.bmp',res=300,width=8,height=8,units='in')

forest(meta6, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "IM more common in MS",
       label.left = "IM more common in Control",
       sortvar = studlab
)
dev.off()


funnel(meta6,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL)

metabias(meta6)

meta10 = update(meta6, byvar = meta_IM$Specific_criteria)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/IM_subgroup_criteria.bmp',res=300,width=8,height=8,units='in')

forest(meta10, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab,
       subgroup = TRUE,
       study.results = TRUE,
)
dev.off()
tiff(file='C:/Users/Ben/Desktop/Work/EBV review/IM_subgroup_criteria_funnel.bmp',res=300,width=8,height=8,units='in')
?tiff
funnel(meta10,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL,
       col=meta_IM$Specific_criteria,
       pch=20,
       cex=1
)
legend(x='topright',c('No specified criteria','Criteria-defined MS'),fill=c('black','red'),cex=1)
dev.off()


#Meta for all seropositive persons (adults and kids)

meta7 = metabin(data=meta_seropositives,
                event.e=meta_seropositives$nMS.,
                n.e=meta_seropositives$nMS,
                event.c=meta_seropositives$nCONTROL.,
                n.c=meta_seropositives$nCONTROL,
                studlab = paste(meta_seropositives$Author,meta_seropositives$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5, 
                subset=(meta_seropositives$serum.CSF=='S')
)

forest(meta7, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)

funnel(meta7,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL)

metabias(meta7)

#criteria-defined MS 

meta9 = update(meta7, byvar=meta_seropositives$technique)
forest(meta9, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab,
       subgroup = TRUE,
       study.results = TRUE
)


meta99 = update(meta7, byvar=meta_seropositives$Adults.paeds)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/seropositivity_subgroup_age.bmp',res=300,width=8,height=8,units='in')

forest(meta99, 
       fontsize = 4, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.4, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab,
       subgroup = TRUE,
       study.results = TRUE,
       bylab=c('Adults','Children'),
       print.byvar = FALSE,
       print.subgroup.labels = TRUE
       
)

dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/seropositives_funnel.bmp',res=300,width=8,height=8,units='in')
funnel(meta99,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL,
       col=meta_seropositives$Adults.paeds,
       pch=20
)
legend(x='topright',c('Adults','Children'),fill=c('black','red'),cex=1)

dev.off()
?funnel

#Meta for PCR detectability
meta11 = metabin(data=meta_PCR,
                event.e=meta_PCR$EBV..MS,
                n.e=meta_PCR$MS.n,
                event.c=meta_PCR$EBV..total.controls,
                n.c=meta_PCR$n.Controls,
                studlab = paste(meta_PCR$Author,meta_PCR$Year,sep=" "),
                method="MH" ,
                sm = 'OR',incr=0.5, 
                subset=(meta_PCR$Sample=='Whole blood' | meta_PCR$Sample=='PBMC')
)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/bloodpcr.bmp',res=300,width=8,height=8,units='in')

forest(meta11, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)
dev.off()
grid::grid.text("Meta-analysis of detectability of EBV by PCR", 0.5, 0.95,gp=grid::gpar(fontsize=14))

funnel(meta11,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL)

metabias(meta11, k.min =5)

#Meta for PCR detectability in CSF
meta12 = metabin(data=meta_PCR,
                 event.e=meta_PCR$EBV..MS,
                 n.e=meta_PCR$MS.n,
                 event.c=meta_PCR$EBV..total.controls,
                 n.c=meta_PCR$n.Controls,
                 studlab = paste(meta_PCR$Author,meta_PCR$Year,sep=" "),
                 method="MH" ,
                 sm = 'OR',incr=0.5, 
                 subset=(meta_PCR$Sample=='CSF' | meta_PCR$Sample== 'CSF (cell-free)' | meta_PCR$Sample=='CSF (cell pellet)')
)

forest(meta12, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)
funnel(meta12,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL)


#Meta for PCR detectability - serum
meta13 = metabin(data=meta_PCR,
                 event.e=meta_PCR$EBV..MS,
                 n.e=meta_PCR$MS.n,
                 event.c=meta_PCR$EBV..total.controls,
                 n.c=meta_PCR$n.Controls,
                 studlab = paste(meta_PCR$Author,meta_PCR$Year,sep=" "),
                 method="MH" ,
                 sm = 'OR',incr=0.5, 
                 subset=(meta_PCR$Sample=='Plasma' | meta_PCR$Sample=='Serum')
)

forest(meta13, 
       fontsize = 6, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.6, 
       squaresize = 0.6,
       col.square = "black",
       label.right = "EBV more common in MS",
       label.left = "EBV more common in Control",
       sortvar = studlab
)
funnel(meta13,
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       level=NULL)


#Import HLA interaction data
meta_HLA_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/meta_HLA_EBV_logodds.csv')
#drop sundquist 2012 due to overlap & drop duplicate simons
meta_HLA_interaction = meta_HLA_interaction[-c(7:10),]


#meta for EBVhiHLA+ OR 
meta14 = metagen(meta_HLA_interaction$lnOR_EBVhiHLApos,
                 meta_HLA_interaction$SElnOR_EBVhiHLApos,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2A.bmp',res=300,width=8,height=8,units='in')
forest(meta14, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta15 = metagen(meta_HLA_interaction$lnOR_EBVloHLApos,
                 meta_HLA_interaction$SElnOR_EBVloHLApos,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2B.bmp',res=300,width=8,height=8,units='in')

forest(meta15, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta16 = metagen(meta_HLA_interaction$lnOR_EBVhiHLAneg,
                 meta_HLA_interaction$SElnOR_EBVhiHLAneg,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2C.bmp',res=300,width=8,height=8,units='in')

forest(meta16, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()

treatment_effects = c(1,exp(meta16$TE.random),exp(meta15$TE.random),exp(meta14$TE.random))
lower_bounds = c(1,exp((meta16$lower.random)),exp((meta15$lower.random)),exp((meta14$lower.random)))
upper_bounds = c(1,exp((meta16$upper.random)),exp((meta15$upper.random)),exp((meta14$upper.random)))
group = factor(c('EBVlo HLA-','EBVhi HLA-','EBVlo HLA+','EBVhi HLA+'),levels = c('EBVlo HLA-','EBVhi HLA-','EBVlo HLA+','EBVhi HLA+'))
ebv_df = data.frame(group,treatment_effects,lower_bounds,upper_bounds)

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2D.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df,mapping=aes(group,treatment_effects, fill=group))+
  geom_col(col='black')+
  theme_classic()+
  labs(y='Odds ratio',x='')+
  geom_errorbar(mapping=aes(group),
                ymin=lower_bounds,
                ymax=upper_bounds,
                width=0.1)+
  scale_y_continuous(limits = c(0,15),breaks=c(0:15))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
                      geom_hline(yintercept=1,lwd=0.5,linetype="dashed")
dev.off()



########################################################
###             3. Interaction analysis              ###
########################################################

### This section plots pooled ORs from HLA-EBNA, HLA-IM, and EBNA-smoking studies

#re Import HLA interaction data
meta_HLA_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/meta_HLA_EBV_logodds.csv')
meta_HLA_interaction = meta_HLA_interaction[meta_HLA_interaction$HLA_method=='Tagging SNPs',]



#meta for EBVhiHLA+ OR 
meta14 = metagen(meta_HLA_interaction$lnOR_EBVhiHLApos,
                 meta_HLA_interaction$SElnOR_EBVhiHLApos,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1A.bmp',res=300,width=8,height=8,units='in')
forest(meta14, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta15 = metagen(meta_HLA_interaction$lnOR_EBVloHLApos,
                 meta_HLA_interaction$SElnOR_EBVloHLApos,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1B.bmp',res=300,width=8,height=8,units='in')

forest(meta15, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta16 = metagen(meta_HLA_interaction$lnOR_EBVhiHLAneg,
                 meta_HLA_interaction$SElnOR_EBVhiHLAneg,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1C.bmp',res=300,width=8,height=8,units='in')

forest(meta16, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()

treatment_effects = c(1,exp(meta16$TE.random),exp(meta15$TE.random),exp(meta14$TE.random))
lower_bounds = c(1,exp((meta16$lower.random)),exp((meta15$lower.random)),exp((meta14$lower.random)))
upper_bounds = c(1,exp((meta16$upper.random)),exp((meta15$upper.random)),exp((meta14$upper.random)))
group = factor(c('EBVlo HLA-','EBVhi HLA-','EBVlo HLA+','EBVhi HLA+'),levels = c('EBVlo HLA-','EBVhi HLA-','EBVlo HLA+','EBVhi HLA+'))
ebv_df = data.frame(group,treatment_effects,lower_bounds,upper_bounds)

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1D.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df,mapping=aes(group,treatment_effects, fill=group))+
  geom_col(col='black')+
  theme_classic()+
  labs(y='Odds ratio',x='')+
  geom_errorbar(mapping=aes(group),
                ymin=lower_bounds,
                ymax=upper_bounds,
                width=0.1)+
  scale_y_continuous(limits = c(0,15),breaks=c(0:15))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
  geom_hline(yintercept=1,lwd=0.5,linetype="dashed")
dev.off()



#re Import HLA interaction data
meta_HLA_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/meta_HLA_EBV_logodds.csv')
meta_HLA_interaction = meta_HLA_interaction[meta_HLA_interaction$HLA_method=='PCR',]
meta_HLA_interaction = meta_HLA_interaction[-4,]


#meta for EBVhiHLA+ OR 
meta14 = metagen(meta_HLA_interaction$lnOR_EBVhiHLApos,
                 meta_HLA_interaction$SElnOR_EBVhiHLApos,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1E.bmp',res=300,width=8,height=8,units='in')
forest(meta14, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta15 = metagen(meta_HLA_interaction$lnOR_EBVloHLApos,
                 meta_HLA_interaction$SElnOR_EBVloHLApos,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1F.bmp',res=300,width=8,height=8,units='in')

forest(meta15, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta16 = metagen(meta_HLA_interaction$lnOR_EBVhiHLAneg,
                 meta_HLA_interaction$SElnOR_EBVhiHLAneg,
                 studlab = paste(meta_HLA_interaction$Author,meta_HLA_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1G.bmp',res=300,width=8,height=8,units='in')

forest(meta16, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()

treatment_effects = c(1,exp(meta16$TE.random),exp(meta15$TE.random),exp(meta14$TE.random))
lower_bounds = c(1,exp((meta16$lower.random)),exp((meta15$lower.random)),exp((meta14$lower.random)))
upper_bounds = c(1,exp((meta16$upper.random)),exp((meta15$upper.random)),exp((meta14$upper.random)))
group = factor(c('EBVlo HLA-','EBVhi HLA-','EBVlo HLA+','EBVhi HLA+'),levels = c('EBVlo HLA-','EBVhi HLA-','EBVlo HLA+','EBVhi HLA+'))
ebv_df = data.frame(group,treatment_effects,lower_bounds,upper_bounds)

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/Supp.fig 1H.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df,mapping=aes(group,treatment_effects, fill=group))+
  geom_col(col='black')+
  theme_classic()+
  labs(y='Odds ratio',x='')+
  geom_errorbar(mapping=aes(group),
                ymin=lower_bounds,
                ymax=upper_bounds,
                width=0.1)+
  scale_y_continuous(limits = c(0,15),breaks=c(0:15))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
  geom_hline(yintercept=1,lwd=0.5,linetype="dashed")
dev.off()


#Import HLA interaction data
meta_HLA_IM_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/IM_HLA_interaction.csv')

#meta for IM and HLA interaction

meta17 = metagen(meta_HLA_IM_interaction$ln.OR._Impos_HLApos,
                 meta_HLA_IM_interaction$SE.lnOR.__Impos_HLApos,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3A.bmp',res=300,width=8,height=8,units='in')

forest(meta17, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab)
dev.off()

meta18 = metagen(meta_HLA_IM_interaction$ln.OR._Impos_HLAneg,
                 meta_HLA_IM_interaction$SE.lnOR._Impos_HLAneg,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3B.bmp',res=300,width=8,height=8,units='in')

forest(meta18, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta19 = metagen(meta_HLA_IM_interaction$ln.OR._Imneg_HLApos,
                 meta_HLA_IM_interaction$SE.lnOR._Imneg_HLApos,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3C.bmp',res=300,width=8,height=8,units='in')
forest(meta19, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()

treatment_effects_2 = c(1,exp(meta17$TE.random),exp(meta18$TE.random),exp(meta19$TE.random))
lower_bounds_2 = c(1,exp((meta17$lower.random)),exp((meta18$lower.random)),exp((meta19$lower.random)))
upper_bounds_2 = c(1,exp((meta17$upper.random)),exp((meta18$upper.random)),exp((meta19$upper.random)))
group_2 = factor(c('IM- HLA-','IM+ HLA+','IM+ HLA-','IM- HLA+'),levels = c('IM- HLA-','IM+ HLA-','IM- HLA+','IM+ HLA+'))
ebv_df_2 = data.frame(group_2,treatment_effects_2,lower_bounds_2,upper_bounds_2)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3D.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df_2,
       mapping=aes(group_2,treatment_effects_2, fill=group_2))+
  geom_col(col='black')+theme_classic()+labs(y='Odds ratio',x='')+
  geom_errorbar(mapping=aes(group_2),
                ymin=lower_bounds_2,
                ymax=upper_bounds_2,
                width=0.1)+
  scale_y_continuous(limits = c(0,15),breaks=c(0:15))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
  geom_hline(yintercept=1,lwd=0.5,linetype="dashed")

dev.off()



#repeat for hla resolution subgroup
#Import HLA interaction data
meta_HLA_IM_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/IM_HLA_interaction.csv')

meta_HLA_IM_interaction = meta_HLA_IM_interaction[c(1,4),]

#meta for IM and HLA interaction

meta17 = metagen(meta_HLA_IM_interaction$ln.OR._Impos_HLApos,
                 meta_HLA_IM_interaction$SE.lnOR.__Impos_HLApos,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/supp.fig 2A.bmp',res=300,width=8,height=8,units='in')

forest(meta17, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab)
dev.off()

meta18 = metagen(meta_HLA_IM_interaction$ln.OR._Impos_HLAneg,
                 meta_HLA_IM_interaction$SE.lnOR._Impos_HLAneg,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/supp.fig2B.bmp',res=300,width=8,height=8,units='in')

forest(meta18, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta19 = metagen(meta_HLA_IM_interaction$ln.OR._Imneg_HLApos,
                 meta_HLA_IM_interaction$SE.lnOR._Imneg_HLApos,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/supp.fig2C.bmp',res=300,width=8,height=8,units='in')
forest(meta19, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()

treatment_effects_2 = c(1,exp(meta17$TE.random),exp(meta18$TE.random),exp(meta19$TE.random))
lower_bounds_2 = c(1,exp((meta17$lower.random)),exp((meta18$lower.random)),exp((meta19$lower.random)))
upper_bounds_2 = c(1,exp((meta17$upper.random)),exp((meta18$upper.random)),exp((meta19$upper.random)))
group_2 = factor(c('IM- HLA-','IM+ HLA+','IM+ HLA-','IM- HLA+'),levels = c('IM- HLA-','IM+ HLA-','IM- HLA+','IM+ HLA+'))
ebv_df_2 = data.frame(group_2,treatment_effects_2,lower_bounds_2,upper_bounds_2)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/suppfig2D.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df_2,
       mapping=aes(group_2,treatment_effects_2, fill=group_2))+
  geom_col(col='black')+theme_classic()+labs(y='Odds ratio',x='')+
  geom_errorbar(mapping=aes(group_2),
                ymin=lower_bounds_2,
                ymax=upper_bounds_2,
                width=0.1)+
  scale_y_continuous(limits = c(0,15),breaks=c(0:15))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
  geom_hline(yintercept=1,lwd=0.5,linetype="dashed")

dev.off()

#repeat for hla resolution subgroup
#Import HLA interaction data
meta_HLA_IM_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/IM_HLA_interaction.csv')

meta_HLA_IM_interaction = meta_HLA_IM_interaction[c(2,3),]

#meta for IM and HLA interaction

meta17 = metagen(meta_HLA_IM_interaction$ln.OR._Impos_HLApos,
                 meta_HLA_IM_interaction$SE.lnOR.__Impos_HLApos,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/supp.fig 2E.bmp',res=300,width=8,height=8,units='in')

forest(meta17, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab)
dev.off()

meta18 = metagen(meta_HLA_IM_interaction$ln.OR._Impos_HLAneg,
                 meta_HLA_IM_interaction$SE.lnOR._Impos_HLAneg,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/supp.fig2F.bmp',res=300,width=8,height=8,units='in')

forest(meta18, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta19 = metagen(meta_HLA_IM_interaction$ln.OR._Imneg_HLApos,
                 meta_HLA_IM_interaction$SE.lnOR._Imneg_HLApos,
                 studlab = paste(meta_HLA_IM_interaction$Author,meta_HLA_IM_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')

bmp(file='C:/Users/Ben/Desktop/Work/EBV review/supp.fig2G.bmp',res=300,width=8,height=8,units='in')
forest(meta19, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()

treatment_effects_2 = c(1,exp(meta17$TE.random),exp(meta18$TE.random),exp(meta19$TE.random))
lower_bounds_2 = c(1,exp((meta17$lower.random)),exp((meta18$lower.random)),exp((meta19$lower.random)))
upper_bounds_2 = c(1,exp((meta17$upper.random)),exp((meta18$upper.random)),exp((meta19$upper.random)))
group_2 = factor(c('IM- HLA-','IM+ HLA+','IM+ HLA-','IM- HLA+'),levels = c('IM- HLA-','IM+ HLA-','IM- HLA+','IM+ HLA+'))
ebv_df_2 = data.frame(group_2,treatment_effects_2,lower_bounds_2,upper_bounds_2)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/suppfig2H.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df_2,
       mapping=aes(group_2,treatment_effects_2, fill=group_2))+
  geom_col(col='black')+theme_classic()+labs(y='Odds ratio',x='')+
  geom_errorbar(mapping=aes(group_2),
                ymin=lower_bounds_2,
                ymax=upper_bounds_2,
                width=0.1)+
  scale_y_continuous(limits = c(0,15),breaks=c(0:15))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
  geom_hline(yintercept=1,lwd=0.5,linetype="dashed")

dev.off()



#Import smoking EBNA interaction data
meta_smoking_interaction = read.csv('C:/Users/Ben/Desktop/Work/EBV review/meta_EBV_smoking_interaction.csv')

#meta for smoking and EBV


meta20 = metagen(meta_smoking_interaction$ln.OR._EBVhiSmokingpos,
                 meta_smoking_interaction$SE.lnOR._EBVhiSmokingpos,
                 studlab = paste(meta_smoking_interaction$Author,meta_smoking_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4A.bmp',res=300,width=8,height=8,units='in')

forest(meta20, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta21 = metagen(meta_smoking_interaction$ln.OR._EBVloSmokingpos,
                 meta_smoking_interaction$SE.lnOR._EBVloSmokingpos,
                 studlab = paste(meta_smoking_interaction$Author,meta_smoking_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4B.bmp',res=300,width=8,height=8,units='in')

forest(meta21, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
meta22 = metagen(meta_smoking_interaction$ln.OR._EBVhiSmokingneg,
                 meta_smoking_interaction$SE.lnOR._EBVhiSmokingneg,
                 studlab = paste(meta_smoking_interaction$Author,meta_smoking_interaction$Year,sep=" "),
                 comb.random = TRUE,
                 sm='OR')
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4C.bmp',res=300,width=8,height=8,units='in')

forest(meta22, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "MS", 
       lab.c = "Control", 
       title="MS", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "MS more likely",
       label.left = "MS less likely",
       leftlabs = c('Study','ln(OR)','SE(ln(OR))','Weight','OR [95% CI]'),
       sortvar = studlab
)
dev.off()
treatment_effects_3 = c(1,exp(meta20$TE.random),exp(meta21$TE.random),exp(meta22$TE.random))
lower_bounds_3 = c(1,exp((meta20$lower.random)),exp((meta21$lower.random)),exp((meta22$lower.random)))
upper_bounds_3 = c(1,exp((meta20$upper.random)),exp((meta21$upper.random)),exp((meta22$upper.random)))
group_3 = factor(c('EBVlo Smoking-','EBVhi Smoking+','EBVlo Smoking+','EBVhi Smoking-'),levels = c('EBVlo Smoking-','EBVlo Smoking+','EBVhi Smoking-','EBVhi Smoking+'))
ebv_df_3 = data.frame(group_3,treatment_effects_3,lower_bounds_3,upper_bounds_3)
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4D.bmp',res=300,width=8,height=8,units='in')
ggplot(data=ebv_df_3,mapping=aes(group_3,treatment_effects_3, fill=group_3))+geom_col(col='black')+theme_classic()+labs(y='Odds ratio',x='')+geom_errorbar(mapping=aes(group_3),ymin=lower_bounds_3,ymax=upper_bounds_3,width=0.1)+
scale_y_continuous(limits = c(0,5),breaks=c(0:5))+
  theme(axis.line = element_line(size=1),
        text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.05,1))+
  scale_fill_manual(values=c(wesanderson::wes_palette('GrandBudapest2',4,type='discrete')))+
  geom_hline(yintercept=1,lwd=0.5,linetype="dashed")
dev.off()




########################################################
###     4. Individual-level Interaction analysis     ###
########################################################


#interaction stuff

#factor 1 - EBV. Factor 2 - HLA. 
HLA_individual = read.csv('C:/Users/Ben/Desktop/Work/EBV review/meta_HLA_EBV_individual.csv')

#factor 1 - EBV. Factor 2 - smoking.
smoking_individual = read.csv('C:/Users/Ben/Desktop/Work/EBV review/smoking_EBV_individual.csv',sep="\t")

#factor 1 - IM. Factor 2 - HLA.
IM_individual = read.csv('C:/Users/Ben/Desktop/Work/EBV review/IM_HLA_individual.csv',sep="\t")


# master function for getting interaction terms of a given table - input is a data.frame of specific format
# this function calculates AP, RERI, S, and multiplicative interaction from raw cross-tabulation
generate_interaction_terms = function(y){
  
#define the populate table function
populate_table = function(x) {
  ebv_table = matrix(ncol=4)
  colnames(ebv_table) = c('Caseness','Riskfactor1','Riskfactor2','Indicator')
  
  for (i in 1:x[[1]]){
    ebv_table = rbind(ebv_table,c(1,0,0,0))
  }
  
  for (i in 1:x[[2]]){
    ebv_table = rbind(ebv_table,c(0,0,0,0))
  }
  
  for (i in 1:x[[3]]){
    ebv_table = rbind(ebv_table,c(1,1,0,1))
  }
  
  for (i in 1:x[[4]]){
    ebv_table = rbind(ebv_table,c(0,1,0,1))
  }
  for (i in 1:x[[5]]){
    ebv_table = rbind(ebv_table,c(1,0,1,2))
  }
  
  for (i in 1:x[[6]]){
    ebv_table = rbind(ebv_table,c(0,0,1,2))
  }
  for (i in 1:x[[7]]){
    ebv_table = rbind(ebv_table,c(1,1,1,3))
  }
  for (i in 1:x[[8]]){
    ebv_table = rbind(ebv_table,c(0,1,1,3))
  }
  ebv_table <<- data.frame(ebv_table)
  
}

#make empty vectors
Attributable_proportion = vector()
SE_Attributable_proportion = vector()
Multiplicative_interaction = vector()
SE_Multiplicative_interaction = vector()
RERI_vector = vector()
SE_RERI_vector = vector()
Synergy_Index = vector()
SE_Synergy_Index = vector()


#define the get_interactions function
get_interaction_terms = function(){
  ebv_table$Indicator = factor(ebv_table$Indicator)
  model_1 = summary(glm(data=data.frame(ebv_table),formula = Caseness ~ Riskfactor1 * Riskfactor2, family = binomial(link="logit")))
  #glm with separately coded interaction vars
  model_2 = summary(glm(data=data.frame(ebv_table),formula = Caseness ~ Indicator, family = binomial(link="logit")))
  model_2
  vcov(model_2)
  covariate_table = vcov(model_2)
  OR_RF1 = exp(model_2$coefficients[2])
  OR_RF2 = exp(model_2$coefficients[3])
  OR_interaction = exp(model_2$coefficients[4])
  RERI = 1 + OR_interaction - OR_RF1 - OR_RF2
  AP <<- RERI / OR_interaction
  Synergy = (OR_interaction-1)/(OR_RF1 + OR_RF2 - 2)
  
  #calculate se of RERI
  h1 = -exp(model_2$coefficients[2])
  h2 = -exp(model_2$coefficients[3])
  h3 = exp(model_2$coefficients[4])
  h1s1 = h1^2*covariate_table[2,2]
  h2s2 = h2^2*covariate_table[3,3]
  h3s3 = h3^2*covariate_table[4,4]
  h1h2s12_2 = 2*h1*h2*covariate_table[2,3]
  h1h3s13_2 = 2*h1*h3*covariate_table[2,4]
  h2h3s23_2 = 2*h2*h3*covariate_table[3,4]
  var_RERI = h1s1 + h2s2 + h3s3 + h1h2s12_2 + h1h3s13_2 + h2h3s23_2
  se_RERI <<- sqrt(var_RERI)
  
        
  # calculate se of AP
  h1 = -exp(model_2$coefficients[2]-model_2$coefficients[4])
  h2 = -exp(model_2$coefficients[3]-model_2$coefficients[4])
  h3 = (exp(model_2$coefficients[2])+exp(model_2$coefficients[3])-1)/exp(model_2$coefficients[4])
  h1s1 = h1^2*covariate_table[2,2]
  h2s2 = h2^2*covariate_table[3,3]
  h3s3 = h3^2*covariate_table[4,4]
  h1h2s12_2 = 2*h1*h2*covariate_table[2,3]
  h1h3s13_2 = 2*h1*h3*covariate_table[2,4]
  h2h3s23_2 = 2*h2*h3*covariate_table[3,4]
  var_AP = h1s1 + h2s2 + h3s3 + h1h2s12_2 + h1h3s13_2 + h2h3s23_2
  se_AP <<- sqrt(var_AP)
  
  #calculate se of synergy
  h1 = -exp(model_2$coefficients[2])/(exp(model_2$coefficients[2])+exp(model_2$coefficients[3])-2)
  h2 = -exp(model_2$coefficients[3])/(exp(model_2$coefficients[2])+exp(model_2$coefficients[3])-2)
  h3 = exp(model_2$coefficients[4])/(exp(model_2$coefficients[4])-1)
  h1s1 = h1^2*covariate_table[2,2]
  h2s2 = h2^2*covariate_table[3,3]
  h3s3 = h3^2*covariate_table[4,4]
  h1h2s12_2 = 2*h1*h2*covariate_table[2,3]
  h1h3s13_2 = 2*h1*h3*covariate_table[2,4]
  h2h3s23_2 = 2*h2*h3*covariate_table[3,4]
  var_synergy = h1s1 + h2s2 + h3s3 + h1h2s12_2 + h1h3s13_2 + h2h3s23_2
  se_synergy <<- sqrt(var_synergy)
  
  
  lnOR_interaction_term = model_1$coefficients[[4]]
  lnOR_interaction_term_se = model_1$coefficients[[8]]
  
  Multiplicative_interaction  <<- c(Multiplicative_interaction,exp(lnOR_interaction_term))
  SE_Multiplicative_interaction <<- c(SE_Multiplicative_interaction,exp(lnOR_interaction_term_se))
  Attributable_proportion <<- c(Attributable_proportion,AP)
  SE_Attributable_proportion <<- c(SE_Attributable_proportion, se_AP)
  RERI_vector <<- c(RERI_vector,RERI)
  SE_RERI_vector <<- c(SE_RERI_vector,se_RERI)
  Synergy_Index <<- c(Synergy_Index,Synergy)
  SE_Synergy_Index <<- c(SE_Synergy_Index,se_synergy)
  
  summary_df <<- data.frame(Multiplicative_interaction,SE_Multiplicative_interaction,Attributable_proportion,SE_Attributable_proportion,RERI_vector,SE_RERI_vector,Synergy_Index,SE_Synergy_Index)
  }

#for loop to fill interaction terms

for (i in 1:nrow(y)){
  populate_table(y[i,c(4:11)])
  get_interaction_terms()
}

study_labels <- c(paste(y$Author,y$Year))

summary_df$log_synergy_index <<- log(Synergy_Index)
meta_AP <<- metagen(summary_df$Attributable_proportion,
                                        summary_df$SE_Attributable_proportion, 
                                        null.effect = 0,
                                        comb.random = TRUE,
                                        studlab = study_labels)

meta_RERI <<- metagen(summary_df$RERI_vector,
                    summary_df$SE_RERI_vector, 
                    null.effect = 0,
                    comb.random = TRUE,
                    studlab = study_labels)


meta_syn <<- metagen(summary_df$log_synergy_index,
                      summary_df$SE_Synergy_Index, 
                      null.effect = log(1),
                      comb.random = TRUE,
                     studlab = study_labels)



meta_Multi_terms <<- metagen(summary_df$Multiplicative_interaction,
                           summary_df$SE_Multiplicative_interaction, 
                                           null.effect = 1,
                                         comb.random = TRUE,
                           studlab = study_labels)



overall_interaction_summary <<- data.frame('Estimate' = c(meta_AP$TE.random,meta_RERI$TE.random,meta_syn$TE.random,meta_Multi_terms$TE.random),
                                           'SE' = c(meta_AP$seTE.random,meta_RERI$seTE.random,meta_syn$seTE.random,meta_Multi_terms$seTE.random),
                                           'P value' = c(meta_AP$pval.random,meta_RERI$pval.random,meta_syn$pval.random,meta_Multi_terms$pval.random)
                                           )

rownames(overall_interaction_summary) <<- c('AP', 'RERI', 'Log(Synergy index)', 'Multiplicative interaction')
}



generate_interaction_terms(HLA_individual[c(1,5),])
write.csv(overall_interaction_summary,file='C:/Users/Ben/Desktop/Work/EBV review/hla_tagging.csv')
generate_interaction_terms(HLA_individual[c(2,3,4),])
write.csv(overall_interaction_summary,file='C:/Users/Ben/Desktop/Work/EBV review/pcr_.csv')
generate_interaction_terms(HLA_individual)
overall_interaction_summary

#generate forest plots
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2H.bmp',res=300,width=8,height=8,units='in')

forest(meta_Multi_terms, 
                        fontsize = 10, 
                        layout ="Revman5", 
                        comb.fixed = FALSE, 
                        lab.e = "Positive interaction", 
                        lab.c = "Negative interaction", 
                        title="AP due to interaction", 
                        spacing = 0.7,
                        ref = 1,
                        xlim='s',
                        at=c(-5,-3,-1,1,3,5,7),
                        squaresize = 0.5,
                        col.square = "black",
                        col.square.lines = 'black',
                        leftlabs = c('Study','Interaction term','SE','Weight','MIT [95% CI]'),
                        label.right = "Positive interaction",
                        label.left = "Negative interaction",
                        sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2G.bmp',res=300,width=8,height=8,units='in')

forest(meta_syn, 
                      fontsize = 10, 
                      layout ="Revman5", 
                      comb.fixed = FALSE, 
                      lab.e = "Positive interaction", 
                      lab.c = "Negative interaction", 
                      title="AP due to interaction", 
                      spacing = 0.7, 
                      ref = log(1),       
                      at=c(-1,0,log(1),1,2,3),
                      xlim=c(-1,log(5)),
                      squaresize = 0.5,
                     leftlabs = c('Study','ln(S)','SE','Weight','ln(S) [95% CI]'),
       
                      col.square = "black",
                      col.square.lines = 'black',
                      label.right = "Positive interaction",
                      label.left = "Negative interaction",
                      sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2F.bmp',res=300,width=8,height=8,units='in')

RERI_forest <<- forest(meta_RERI, 
                       fontsize = 10, 
                       layout ="Revman5", 
                       comb.fixed = FALSE, 
                       lab.e = "Positive interaction", 
                       lab.c = "Negative interaction", 
                       title="AP due to interaction", 
                       spacing = 0.7, 
                       squaresize = 0.5,
                       xlim=c(-10,10),
                       leftlabs = c('Study','RERI','SE','Weight','RERI [95% CI]'),
                       
                       col.square = "black",
                       col.square.lines = 'black',
                       label.right = "Positive interaction",
                       label.left = "Negative interaction",
                       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_2E.bmp',res=300,width=8,height=8,units='in')

forest(meta_AP, 
                     fontsize = 10, 
                     layout ="Revman5", 
                     comb.fixed = FALSE, 
                     lab.e = "Positive interaction", 
                     lab.c = "Negative interaction", 
                     title="AP due to interaction", 
                     spacing = 0.7, 
                     squaresize = 0.5,
                     col.square = "black",
                     col.square.lines = 'black',
                     label.right = "Positive interaction",
                     label.left = "Negative interaction",
                     sortvar = studlab,
                      leftlabs = c('Study','AP','SE','Weight','AP [95% CI]'),
       
)
dev.off()



overall_interaction_summary



generate_interaction_terms(IM_individual)
#generate forest plots
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3H.bmp',res=300,width=8,height=8,units='in')

forest(meta_Multi_terms, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "Positive interaction", 
       lab.c = "Negative interaction", 
       title="AP due to interaction", 
       spacing = 0.7,
       ref = 1,
       xlim='s',
       at=c(-5,-3,-1,1,3,5,7),
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       leftlabs = c('Study','Interaction term','SE','Weight','MIT [95% CI]'),
       label.right = "Positive interaction",
       label.left = "Negative interaction",
       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3G.bmp',res=300,width=8,height=8,units='in')

forest(meta_syn, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "Positive interaction", 
       lab.c = "Negative interaction", 
       title="AP due to interaction", 
       spacing = 0.7, 
       ref = log(1),       
       at=c(-1,0,log(1),1,2,3),
       xlim=c(-1,log(5)),
       squaresize = 0.5,
       leftlabs = c('Study','ln(S)','SE','Weight','ln(S) [95% CI]'),
       col.square = "black",
       col.square.lines = 'black',
       label.right = "Positive interaction",
       label.left = "Negative interaction",
       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3F.bmp',res=300,width=8,height=8,units='in')

RERI_forest <<- forest(meta_RERI, 
                       fontsize = 10, 
                       layout ="Revman5", 
                       comb.fixed = FALSE, 
                       lab.e = "Positive interaction", 
                       lab.c = "Negative interaction", 
                       title="AP due to interaction", 
                       spacing = 0.7, 
                       squaresize = 0.5,
                       xlim=c(-10,10),
                       leftlabs = c('Study','RERI','SE','Weight','RERI [95% CI]'),
                       
                       col.square = "black",
                       col.square.lines = 'black',
                       label.right = "Positive interaction",
                       label.left = "Negative interaction",
                       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_3E.bmp',res=300,width=8,height=8,units='in')

forest(meta_AP, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "Positive interaction", 
       lab.c = "Negative interaction", 
       title="AP due to interaction", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "Positive interaction",
       label.left = "Negative interaction",
       sortvar = studlab,
       leftlabs = c('Study','AP','SE','Weight','AP [95% CI]'),
       
)
dev.off()

generate_interaction_terms(IM_individual[-1,])
overall_interaction_summary




generate_interaction_terms(smoking_individual)
overall_interaction_summary
#generate forest plots
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4H.bmp',res=300,width=8,height=8,units='in')

forest(meta_Multi_terms, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "Positive interaction", 
       lab.c = "Negative interaction", 
       title="AP due to interaction", 
       spacing = 0.7,
       ref = 1,
       xlim='s',
       at=c(-5,-3,-1,1,3,5,7),
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       leftlabs = c('Study','Interaction term','SE','Weight','MIT [95% CI]'),
       label.right = "Positive interaction",
       label.left = "Negative interaction",
       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4G.bmp',res=300,width=8,height=8,units='in')

forest(meta_syn, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "Positive interaction", 
       lab.c = "Negative interaction", 
       title="AP due to interaction", 
       spacing = 0.7, 
       ref = log(1),       
       at=c(-1,0,log(1),1,2,3),
       xlim=c(-1,log(5)),
       squaresize = 0.5,
       leftlabs = c('Study','ln(S)','SE','Weight','ln(S) [95% CI]'),
       col.square = "black",
       col.square.lines = 'black',
       label.right = "Positive interaction",
       label.left = "Negative interaction",
       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4F.bmp',res=300,width=8,height=8,units='in')

RERI_forest <<- forest(meta_RERI, 
                       fontsize = 10, 
                       layout ="Revman5", 
                       comb.fixed = FALSE, 
                       lab.e = "Positive interaction", 
                       lab.c = "Negative interaction", 
                       title="AP due to interaction", 
                       spacing = 0.7, 
                       squaresize = 0.5,
                       xlim=c(-10,10),
                       leftlabs = c('Study','RERI','SE','Weight','RERI [95% CI]'),
                       
                       col.square = "black",
                       col.square.lines = 'black',
                       label.right = "Positive interaction",
                       label.left = "Negative interaction",
                       sortvar = studlab
)
dev.off()
bmp(file='C:/Users/Ben/Desktop/Work/EBV review/figure_4E.bmp',res=300,width=8,height=8,units='in')

forest(meta_AP, 
       fontsize = 10, 
       layout ="Revman5", 
       comb.fixed = FALSE, 
       lab.e = "Positive interaction", 
       lab.c = "Negative interaction", 
       title="AP due to interaction", 
       spacing = 0.7, 
       squaresize = 0.5,
       col.square = "black",
       col.square.lines = 'black',
       label.right = "Positive interaction",
       label.left = "Negative interaction",
       sortvar = studlab,
       leftlabs = c('Study','AP','SE','Weight','AP [95% CI]'),
       
)
dev.off()

generate_interaction_terms(smoking_individual)
generate_interaction_terms(smoking_individual[-1,])
write.csv(overall_interaction_summary,file='C:/Users/Ben/Desktop/Work/EBV review/smoking_nosecondhand.csv')
generate_interaction_terms(smoking_individual[-2,])
write.csv(overall_interaction_summary,file='C:/Users/Ben/Desktop/Work/EBV review/smoking_nocot.csv')


