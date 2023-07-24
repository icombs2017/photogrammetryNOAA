if (!require("pacman")) install.packages("pacman")
pacman::p_load("ggplot2","officer","ggpubr", "rcompanion", "RColorBrewer", 
               "patchwork", "magrittr","reshape2", "stringr", "plyr", "dplyr", 
               "flextable", "tidyr", "tibble", "vegan", "readxl")


acer <- readxl::read_xlsx("../data/TagLab/ACERpercentCover.xlsx", sheet = 1, col_names = TRUE, na = c("", "N/A", "NA", "NA NA"))
acer$site <- as.factor(acer$site)
acer$reef <- as.factor(acer$reef)
acer$date <- as.factor((acer$date))




acerPlot1 <- ggplot(acer, aes(x = date, y = percentCoverAcer, fill = site))+
  geom_bar(stat = 'identity', color = 'black',)+
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(breaks = seq(0,100,by = 20))+
  scale_x_discrete(labels = c("T0", "T12"))+
  ggtitle(expression(paste("Percent Coverage of ", italic("Acropora cervicornis"))))+
  labs(x = "Time", fill = "Site", y = "Percent Coverage (%)")+
  facet_wrap(~reef, scales = "free")

acerPlot <- acerPlot1+theme(
  axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5, hjust = 0.5, face= "bold"), 
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 40, face = "bold"), 
  axis.text.y = element_text(colour = "black", size = 40, face = "bold"),
  legend.title = element_text(size = 40, face = "bold"), 
  legend.text = element_text(size = 30, colour = "black"), 
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
  panel.background = element_rect(fill = '#F5F5F5'), 
  plot.title = element_text(size = 40, face = "bold"), 
  axis.line = element_line(colour = "black"), 
  axis.ticks = element_line(color="black"), 
  text = element_text(size=40, color="black"), 
  legend.position = "right")
acerPlot
  
ggsave("../figures/percentCoverAcer.png", plot = acerPlot, width = 20, height = 15, units = "in", dpi = 600)


