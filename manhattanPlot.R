

library("optparse")
# ==========================================================================

# ==========================================================================
option_list = list(
  make_option(c("--trait"), type = "character"),
  make_option(c("--omit"), type = "character"),
  make_option(c("--anc"), type = "character"),  
  make_option(c("--out"), type = "character")  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)


library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)

options(datatable.fread.datatable=FALSE)


sumstat = fread(paste0("metal_geshtn_European.African.Asian.Hispanic_allBiobanks_omitNone_230203_1.txt")) # with columns Chromosome, Position, P, Freq1, P-value


sumstat = sumstat %>%
	rename(c("P"="P-value"))

sumstat$Freq1 = as.numeric(sumstat$Freq1)
sumstat$P = as.numeric(sumstat$P)
sumstat$Chromosome = as.numeric(sumstat$Chromosome)
sumstat$Position = as.numeric(sumstat$Position)

#################################


sumstat = sumstat %>%
	filter(!(is.na(Freq1) | is.na(P) | is.na(Chromosome) | is.na(Position)))
dim(sumstat)


ff = 0.01

sumstat_common = sumstat %>%
	filter(Freq1 >= ff & Freq1 <= 1-ff) %>%
	filter(Chromosome <= 22 & Chromosome>=1)
	

dim(sumstat_common)

ff1 = 0.4
ff2 = 0.01
sumstat_common_p = sumstat_common %>% 
	filter((P <= 5e-8 & (Freq1 >= ff2 & Freq1 <= 1-ff2)) | (Freq1 >= ff1 & Freq1 <= 1-ff1))

a = sumstat_common_p %>% filter(Chromosome==3)
dim(a)

dim(sumstat_common_p)

# sumstat_common_p$Position = as.numeric(sumstat_common_p$Position)


sumstat_common_s = sumstat_common_p %>%
	group_by(Chromosome) %>%
	summarise(chr_len=max(Position)) %>%
	mutate(tot=cumsum(chr_len)-chr_len) %>%
	left_join(sumstat_common_p, ., by=c("Chromosome"="Chromosome")) %>%
	arrange(Chromosome, Position) %>%
	mutate(Positioncum=Position+tot)

axisdf = sumstat_common_s %>% 
	group_by(Chromosome) %>% 
	summarize(center=( max(Positioncum) + min(Positioncum) ) / 2 )
axisdf


sumstat_common_s %>% 
		group_by(Chromosome) %>%
		summarise(chr_min=min(Position), chr_max=max(Position))



sumstat_common_s1 = sumstat_common_s
sumstat_common_s1$color = NA
sumstat_common_s1$gene = NA


pos_plot = fread("pos_plot_geshtn.txt") ## if user has candidate lead SNPs/genes to plot with columns chr, pos, Gene

# pos_plot = pos_plot %>% filter(Gene != "MECOM")
pos_plot = pos_plot %>% filter(chr!=3)


for (i in 1:nrow(pos_plot)) {
	row = pos_plot[i,]
	rsid = paste0(row$chr, ":", row$pos)
	
	idx = which(sumstat_common_s1$MarkerName == rsid)
    print(row)
    print(idx)
    if (row$Gene != "dd") sumstat_common_s1$gene[idx] = row$Gene
    
	idx1 = which(sumstat_common_s1$Chromosome == row$chr & sumstat_common_s1$Position >= row$pos-500000 & sumstat_common_s1$Position <= row$pos+500000)
	sumstat_common_s1$color[idx1] = "red"
}

sumstat_common_s1 = sumstat_common_s1 %>%
    filter(!(P<=5e-8 & color!="red"))

sumstat_common_s_red = sumstat_common_s1 %>% filter(color=="red")
sumstat_common_s_red = sumstat_common_s_red[order(sumstat_common_s_red$P),]
sumstat_common_sp = sumstat_common_s1

dim(sumstat_common_sp)



a = sumstat_common_s %>% filter(Chromosome==3)
dim(a)
a[which.min(a$P),]

reso <- 300
length <- 3.25*reso/72
min(sumstat_common_s$P)

png(paste0("AllBiobanks_geshtn_230502.png"),units="in", res=reso,height=5,width=15, pointsize=2)

ggplot(sumstat_common_sp, aes(x=Positioncum, y=-log10(P))) +
	geom_point( aes(color=as.factor(Chromosome)), alpha=0.8, size=1.3) +
	scale_color_manual(values = rep(c("#AAAAAA", "#BBBBBB"), 22 )) +
	# ylim(0, 17.5) +
	geom_point(data=sumstat_common_s_red, color="red2", alpha=0.8, size=1.3) +
    scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
	scale_y_continuous(expand = c(0, 0), limits=c(0, 17.5)) +     # remove space between plot area and x axis
	geom_hline(yintercept=-log10(5e-8), color="brown") +
	# Custom the theme:
	theme_bw(base_size=15) +
	theme( 
	    legend.position="none",
	    panel.border = element_blank(),
	    panel.grid.major.x = element_blank(),
	    panel.grid.minor.x = element_blank()
	) +
	labs(x="")

dev.off()


