library(ggplot2)
library(patchwork)
library(stats)
library(viridis)
library(ggpointdensity)
library(ggrastr)


setwd("~/yeast-idr-analysis-main/processed_scores/")
wt = read.csv("wt_scores_sortSeq.csv")
wt = wt[wt$iRFP_mean >= -1,] #select for only stable peps 
wt = wt[,c(2:4)] #remove stability scores for downstream analysis


ccr4 = read.csv("ccr4_scores.csv")
pop2 = read.csv('pop2_scores.csv')
dcp2 = read.csv('dcp2_scores.csv')


fit_loess <-function(df_1=wt,df_2, cond){
  #predict loess
  combo = merge(df_1, df_2, by.x = 'Peptide', by.y='Name')
  mod_1 = loess(combo$YFP_mean~combo$activity_score_r1, span=0.25, control=loess.control(surface="direct"))
  mod_2 = loess(combo$YFP_mean~combo$activity_score_r2, span=0.25, control=loess.control(surface="direct"))

  #map loess activity and get mean/sd
  combo$loess_rep1 = predict(mod_1)
  combo$loess_rep2 = predict(mod_2)
  combo$avg_loess = rowMeans(combo[,c("loess_rep1",'loess_rep2')])
  combo$loess_sd = apply(combo[,c("loess_rep1",'loess_rep2')],1,sd)
  
  #compute Z-test
  global_sd = sqrt(mean(combo$YFP_sd)^2 + mean(combo$loess_sd)^2)/sqrt(2)
  combo$diff = combo$avg_loess - combo$YFP_mean
  combo$z_val_raw = sapply(seq(1,nrow(combo)), function(x){2*pnorm(-abs(combo$diff[x])/global_sd)})
  combo$z_adj = p.adjust(combo$z_val_raw, 'BH')
  
  #now get magnitude corresponding to sig diff
  targ_z_adj = max(combo$z_adj[combo$z_adj < 0.05]) #find max z_adj that still is sign
  targ_z_raw = combo[which.min(abs(combo$z_adj - targ_z_adj)),'z_val_raw'] #get the raw pval this correpsonds
  mag = qnorm(targ_z_raw/2)* global_sd #get the difference magnitude this corresponds to 
  
  #now return new column that predicts the loess lines demarcating sig
  combo$xes = seq(-2,2,length.out = nrow(combo))
  combo$y1 = predict(mod_1, newdata=combo$xes) ###there is something buggy here with mod_1 and prediction; will figure out later
  
  
  combo$y2 = predict(mod_2, newdata=combo$xes)
  combo$y_mean_above = rowMeans(combo[,c('y1','y2')])+mag
  combo$y_mean_below = rowMeans(combo[,c('y1','y2')])-mag
  colnames(combo)[4:ncol(combo)] = paste0(cond,'_',colnames(combo)[4:ncol(combo)]) #format colnames
  
  return(combo)
  
}

#fit all to the loess and get padj
ccr4_fit = fit_loess(df_2= ccr4, cond='ccr4')
pop2_fit = fit_loess(df_2=pop2, cond='pop2')
dcp2_fit = fit_loess(df_2=dcp2, cond='dcp2')

#create density plots
ccr4_plot = ggplot(ccr4_fit, aes(x=YFP_mean,y=ccr4_avg_activity))+
  geom_pointdensity(size = 0.5, show.legend = FALSE) + scale_color_viridis() + 
  geom_line(aes(x=ccr4_y_mean_above,y=ccr4_xes),color='#bdbdbd') +
  geom_line(aes(x=ccr4_y_mean_below,y=ccr4_xes),color='#bdbdbd') +
  scale_y_continuous( expand=c(0,0), limits = c(-2,2)) + 
  scale_x_continuous(expand=c(0,0), limits = c(-2,2)) + coord_fixed() + 
  theme_classic() + ylab('Ccr4 Activity Score') + xlab('WT Activity Score') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


pop2_plot = ggplot(pop2_fit, aes(x=YFP_mean,y=pop2_avg_activity))+
  geom_pointdensity(size =0.5, show.legend = FALSE) + scale_color_viridis() + 
  geom_line(aes(x=pop2_y_mean_above,y=pop2_xes),color='#bdbdbd') +
  geom_line(aes(x=pop2_y_mean_below,y=pop2_xes),color='#bdbdbd') +
  scale_y_continuous( expand=c(0,0), limits = c(-2,2)) + 
  scale_x_continuous(expand=c(0,0), limits = c(-2,2)) + coord_fixed() + 
  theme_classic() + ylab('Pop2 Activity Score') + xlab('WT Activity Score') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

dcp2_plot = ggplot(dcp2_fit, aes(x=YFP_mean,y=dcp2_avg_activity))+
  geom_pointdensity(size=0.5, show.legend = FALSE) + scale_color_viridis() + 
  geom_line(aes(x=dcp2_y_mean_above,y=dcp2_xes),color='#bdbdbd') +
  geom_line(aes(x=dcp2_y_mean_below,y=dcp2_xes),color='#bdbdbd') +
  scale_y_continuous( expand=c(0,0), limits = c(-2,2)) + 
  scale_x_continuous(expand=c(0,0), limits = c(-2,2)) + coord_fixed() + 
  theme_classic() + ylab('Dcp2 Activity Score') + xlab('WT Activity Score') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


ccr4_plot = rasterize(ccr4_plot, layers='Point',dpi=300)
pop2_plot = rasterize(pop2_plot, layers='Point',dpi=300)
dcp2_plot = rasterize(dcp2_plot, layers='Point',dpi=300)

composite_plot <- ccr4_plot + pop2_plot + dcp2_plot
composite_plot <- composite_plot + plot_layout(ncol = 3)

setwd("~/yeast-idr-analysis-main/fig3_ms1/")
if(!file.exists('wt_vs_maid_actiivty.pdf')){
  ggsave('wt_vs_maid_actiivty.pdf', composite_plot,height=5, width=20, units='cm')
}


###get number of repressors that increase activity
#145 for Ccr4, 331 for Pop2, 277 for Dcp2
ccr4_reps = ccr4_fit[ccr4_fit$YFP_mean <= -1 & ccr4_fit$ccr4_diff >= 0 & ccr4_fit$ccr4_z_adj < 0.05,]
pop2_reps = pop2_fit[pop2_fit$YFP_mean <= -1 & pop2_fit$pop2_diff >= 0 & pop2_fit$pop2_z_adj < 0.05,]
dcp2_reps = dcp2_fit[dcp2_fit$YFP_mean <= -1 & dcp2_fit$dcp2_diff >= 0 & dcp2_fit$dcp2_z_adj < 0.05,]


##now, find frags that change activity between conditions
library(UpSetR)
#merge by peptide name and get the z_adj 
master <- merge(merge(ccr4_fit, pop2_fit[c(1,4:ncol(pop2_fit))], by = "Peptide", all = FALSE), 
                dcp2_fit[c(1,4:ncol(dcp2_fit))], by = "Peptide", all = FALSE)
#select repressors only
master = master[master$YFP_mean <= -1, ]

upset_input = data.frame('Peptide'=master$Peptide,
                         'ccr4' = ifelse(master$ccr4_z_adj < 0.05 & master$ccr4_diff > 0,1,0),
                         'pop2' = ifelse(master$pop2_z_adj < 0.05 & master$pop2_diff > 0,1,0),
                         'dcp2' = ifelse(master$dcp2_z_adj < 0.05 & master$dcp2_diff > 0,1,0))


#make upset plot, manually saved
upset(upset_input, order.by = 'freq')

#number of genes that don't change in any background
num_diff = length(which(rowSums(upset_input[,2:4]) > 0))
num_noLoss = length(which(rowSums(upset_input[,2:4]) == 0))

######plot indiv protein comparison across conditions
setwd("~/yeast-idr-analysis-main/processed_scores/")
protLen = read.csv("protein_len.csv")

master = merge(merge(ccr4_fit, pop2_fit, by='Peptide'),dcp2_fit, by='Peptide')
master$yorf = sapply(master$Peptide, function(x){sub("^(.*?)_.*$", "\\1", x)})
master$start = as.numeric(sapply(master$Peptide, function(x){sub("^.*?_(.*?)_.*$", "\\1", x)}))
master$stop = as.numeric(sapply(master$Peptide, function(x){sub("^.*_(.*)$", "\\1", x)}))
master$len = as.numeric(sapply(master$yorf, function(x){protLen[match(x, protLen$Standard.Name),'Length']}))

plot_epi = function(df=master,gene_name, spacing=50){
  gene_df = df[df$yorf==gene_name,]
  gene_df = gene_df[order(gene_df$start),]
  gene_df$group = cumsum(c(gene_df$start[1], diff(gene_df$start) > spacing))
  return(gene_df)
}

edc1 = plot_epi(gene_name = 'Edc1p')
edc1_zoom = edc1[c(6,8),]

cols = c('#66c2a5','#fc8d62','#8da0cb')

mot_1 = data.frame(xmin=153, xmax=158, ymin=-Inf,ymax=Inf)
mot_2 = data.frame(xmin=169, xmax=172, ymin=-Inf,ymax=Inf)


edc1_epi <- ggplot(edc1_zoom, aes(x = start, y = YFP_mean.x)) +
  geom_rect(data=mot_1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#bdbdbd", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=mot_2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#bdbdbd", alpha=0.5, inherit.aes = FALSE) +
  geom_segment(aes(x=start,xend = stop, y= YFP_mean.x, yend = YFP_mean.x), color='black', size=2) + 
  geom_segment(aes(x=start,xend = stop, y= ccr4_avg_loess, yend = ccr4_avg_loess),color=cols[1], size=2) + 
  geom_segment(aes(x=start,xend = stop, y= pop2_avg_loess, yend = pop2_avg_loess),color=cols[2], size=2) + 
  geom_segment(aes(x=start,xend = stop, y= dcp2_avg_loess, yend = dcp2_avg_loess),color=cols[3], size=2) +
  geom_errorbar(aes(x = start+24, ymin = YFP_mean.x - YFP_sd.x, ymax = YFP_mean.x + YFP_sd.x), width = 3, color='black',size=0.3) +
  geom_errorbar(aes(x = start+24, ymin = ccr4_avg_loess - ccr4_loess_sd,
                    ymax = ccr4_avg_loess + ccr4_loess_sd), width = 3, color='black',size=0.3) +
  geom_errorbar(aes(x = start+24, ymin = pop2_avg_loess - pop2_loess_sd,
                    ymax = pop2_avg_loess + pop2_loess_sd), width = 3, color='black',size=0.3) +
  geom_errorbar(aes(x = start+24, ymin = dcp2_avg_loess - dcp2_loess_sd,
                    ymax = dcp2_avg_loess + dcp2_loess_sd), width = 3, color='black',size=0.3) +
  ylim(-2,1) + 
  ylab('Activity Score') + xlab('Residue') + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/yeast-idr-analysis-main/fig3_ms1/edc1_zoom/")
if(!file.exists('edc1_zoom.pdf')){
  ggsave('edc1_zoom.pdf',edc1_epi,width = 8, height=5, units='cm')
}
