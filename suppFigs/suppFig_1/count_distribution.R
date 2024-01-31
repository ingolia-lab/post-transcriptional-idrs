library(ggplot2)
setwd("~/yeast-idr-analysis-main/processed_scores/")
scores = read.csv("wt_yfp_iRFP_counts.csv")

#format data
pho = scores[scores$Name == 'Pho92p_17_66',]
pho_1 = data.frame(bin=seq(1,4), reads=t(pho[,c('yfp_FL_1','yfp_L_1','yfp_R_1','yfp_FR_1')]))
pho_1[,2] = pho_1[,2]/sum(pho_1[,2])
pho_2 = data.frame(bin=seq(1,4), reads=t(pho[,c('yfp_FL_2','yfp_L_2','yfp_R_2','yfp_FR_2')]))
pho_2[,2] = pho_2[,2]/sum(pho_2[,2])
colnames(pho_1) = c('bins', 'norm_reads')
colnames(pho_2) = c('bins', 'norm_reads')


gcd = scores[scores$Name == 'Gcd11p_11_60',]
gcd_1 = data.frame(bin=seq(1,4), reads=t(gcd[,c('yfp_FL_1','yfp_L_1','yfp_R_1','yfp_FR_1')]))
gcd_1[,2] = gcd_1[,2]/sum(gcd_1[,2])
gcd_2 = data.frame(bin=seq(1,4), reads=t(gcd[,c('yfp_FL_2','yfp_L_2','yfp_R_2','yfp_FR_2')]))
gcd_2[,2] = gcd_2[,2]/sum(gcd_2[,2])
colnames(gcd_1) = c('bins', 'norm_reads')
colnames(gcd_2) = c('bins', 'norm_reads')

pal = c('#ca0020','#f4a582','#92c5de','#0571b0')

read_dist <- ggplot() + 
  geom_point(data = pho_1, aes(x = bins, y = norm_reads), size = 2, color = pal[1]) + 
  geom_line(data = pho_1, aes(x = bins, y = norm_reads), color = pal[1]) + 
  geom_point(data = pho_2, aes(x = bins, y = norm_reads), size = 2, color = pal[2]) + 
  geom_line(data = pho_2, aes(x = bins, y = norm_reads), color = pal[2]) + 
  geom_point(data = gcd_1, aes(x = bins, y = norm_reads), size = 2, color = pal[3]) + 
  geom_line(data = gcd_1, aes(x = bins, y = norm_reads), color = pal[3]) + 
  geom_point(data = gcd_2, aes(x = bins, y = norm_reads), size = 2, color = pal[4]) + 
  geom_line(data = gcd_2, aes(x = bins, y = norm_reads), color = pal[4]) + 
  theme_classic() + ylab('Normalized Reads') + xlab('Bin') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/yeast-idr-analysis-main/suppFigs/suppFig_1/")
ggsave('pho_gcd_counts.pdf',read_dist, height = 5, width=5, units='cm' )

