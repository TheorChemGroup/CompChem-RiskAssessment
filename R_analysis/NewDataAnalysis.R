library(tidyverse)
library(ggrepel)
library(pdqr)
library(reshape2)
library(HDInterval)
library(ggridges)

ForthPropBayes <- function(Product,Energy,SDdistr,temp,runs=10000,lowst=F,reaction="something"){
  Temp <- unique(temp)
  print(paste("Working on",unique(reaction)))
  if(lowst){
    InT <- data_frame(Product,Energy) %>%
      group_by(Product) %>% 
      summarise(Energy=min(Energy,na.rm=T)) %>% 
      crossing(run = 1:runs) %>% 
      mutate(noise = rnorm(n(), sd = SDdistr(n())))
  }
  else{
    InT <- data_frame(Product,Energy) %>% 
      crossing(run = 1:runs) %>% 
      mutate(noise = rnorm(n(), sd = SDdistr(n())))
  }
  time <- system.time(InT2 <- InT %>% 
                        group_by(run) %>% 
                        summarise(ee_theor_noisy=NRJ2ee2(Energy + noise,Product,Temp,"R",lowst)) %>% 
                        .$ee_theor_noisy
  )
  print(time[3])
  print(InT2)
  return(InT2)
}

NRJ2ee2 <- function(Energy,Product,temp, relat = "R",lowest=F){
  Temp <- unique(temp)
  if(lowest){
    TSc=exp(-1000*(min(Energy[Product!=relat],na.rm=T)-min(Energy[Product==relat],na.rm=T))/(1.9872036*Temp))
    100*(1-TSc)/(1+TSc)
  }
  else {
    TSc=exp((-1000*(Energy))/(1.9872036*Temp))
    TSc=100*TSc/sum(TSc,na.rm=T)
    sum(TSc[Product==relat],na.rm=T) * 2 - 100
  }
}

AllData <- read_delim("AllData.csv",";")

AllData %>% filter(Cat=="All") %>%
  #filter(Solv=="PCM", DFT=="B97D", EnType=="E") %>% 
  group_by(Solv, DFT, EnType) %>% 
  summarize(sigma_mlh=SD[which.min(abs(SDLH-1))])

Figure2B <- AllData %>% 
  filter(!Cat%in%c("1a-4e","All","AllRuined")) %>% 
  ggplot(aes(x=SD,color=Cat,y=SDLH)) +
  geom_line() +
  geom_line(data=AllData %>% filter(Cat=="All", Solv=="PCM", DFT=="B97D", EnType=="E"), 
            size=0.7, color='#e41a1c', alpha=0.8) +
  geom_line(data=AllData %>% filter(Cat=="1a-4e", Solv=="PCM", DFT=="B97D", EnType=="E"), 
            size=0.7, color='#e41a1c', alpha=0.8, linetype="dotted") +
  theme_bw() +
  coord_cartesian(expand = 0, xlim=c(0,1.25), ylim=c(0,1.05)) +
  scale_y_continuous(name="Likelihood",breaks=NULL) +
  scale_x_continuous(name="σ, kcal/mol",limits=c(0,1.01),breaks=c(0.2,0.4,0.6,0.8,1.0)) +
  geom_text_repel(data=AllData %>% 
                    filter(!Cat%in%c("1a-4e","All","AllRuined")) %>%
                    group_by(Cat) %>%
                    filter(SD==SD[which.min(abs(SD-1.0))]), 
                  aes(label=Cat,x=1,y=SDLH), fontface="bold",
                  xlim=c(1,1.25), arrow=arrow(angle=15,length=unit(0.5,"mm")), min.segment.length=0,
                  max.overlaps = Inf, size=3, force_pull = 20, force = 0.1, box.padding=0.2, segment.alpha=0.5, segment.size=0.3) +
  theme(legend.position="none",
        axis.title.y = element_text(margin=margin(0,0,0,0,"mm")),
        panel.grid = element_blank())
Figure2B  
ggsave("Figure2B.png",Figure2B,width=8,height=5,dpi=2000, units="cm")

Figure2C <- AllData %>% 
  filter(Cat%in%c("All"),EnType=="E") %>% 
  ggplot(aes(x=SD,y=SDLH,color=DFT,linetype=Solv)) +
  scale_color_brewer(palette = "Set1",name="DFT") +
  scale_linetype_manual(limits=c("PCM","SMD","None"),
                        name="Solv.\nmodel",
                        values = c("solid","dashed","dotted")) +
  scale_y_continuous(name="Likelihood",limits=c(0,1.1),breaks=c(0.0,1.0)) +
  scale_x_continuous(name="σ, kcal/mol",limits=c(0.15,1.0),breaks=c(0.2,0.4,0.6,0.8,1.0)) +
  geom_line() +
  theme_bw() +
  theme(legend.spacing = unit(0,units="mm"),
        legend.margin = margin(2,0,0,0,"mm"),
        legend.box.margin = margin(5,0,0,0,"mm"),
        legend.key.height = unit(3,"mm"),
        axis.title.y = element_text(margin=margin(0,0,0,0,"mm"),vjust=-3),
        panel.grid = element_blank())
Figure2C  
ggsave("Figure2C.png",Figure2C,width=8,height=4,dpi=2000, units="cm")

FigureS1 <- AllData %>% 
  filter(Cat%in%c("All"),DFT=="B97D",Solv=="PCM") %>% 
  ggplot(aes(x=SD,y=SDLH,color=EnType)) +
  scale_color_brewer(palette = "Set1",name="Energy type") +
  scale_y_continuous(name="Likelihood",limits=c(0,1.1),breaks=c(0.0,1.0)) +
  scale_x_continuous(name="σ, kcal/mol",limits=c(0,1.0),breaks=c(0.2,0.4,0.6,0.8,1.0)) +
  geom_line() +
  theme_bw() +
  theme(legend.spacing = unit(1,units="mm"),
        legend.margin = margin(0,0,0,0,"mm"),
        legend.key.height = unit(4.5,"mm"),
        axis.title.y = element_text(margin=margin(0,0,0,0,"mm"),vjust=-3),
        panel.grid = element_blank())
FigureS1
ggsave("FigureS1.png",FigureS1,width=8,height=6,dpi=2000, units="cm")

Figure2C2 <- AllData %>% 
  filter(Cat%in%c("All","AllRuined"),EnType=="E", Solv%in%c("None","PCM")) %>%
  mutate(Method=ifelse(Solv=="None",DFT,paste0(DFT,"/",Solv))) %>% 
  ggplot(aes(x=SD,y=SDLH,color=Method,linetype=Cat)) +
  scale_color_brewer(palette = "Set2",name="Method") +
  scale_linetype_manual(limits=c("AllRuined","All"),
                        name="Missed\ntransition\nstate",
                        labels=c("Yes","No"),
                        values = c("solid","dashed")) +
  scale_y_continuous(name="Likelihood",limits=c(0,1.1),breaks=c(0.0,1.0)) +
  scale_x_continuous(name="σ, kcal/mol",limits=c(0.15,3.0),breaks=seq(0,3,0.5)) +
  geom_line() +
  theme_bw() +
  theme(legend.spacing = unit(0,units="mm"),
        legend.margin = margin(2,0,0,0,"mm"),
        legend.box.margin = margin(5.5,0,0,0,"mm"),
        legend.key.height = unit(3,"mm"),
        axis.title.y = element_text(margin=margin(0,0,0,0,"mm"),vjust=-3),
        panel.grid = element_blank())
Figure2C2  
ggsave("Figure2C2.png",Figure2C2,width=9.4,height=3.84,dpi=2000, units="cm")

##### Preparing distributions for comparison ----

SDDistrTab <- AllData %>% filter(Cat=="All") %>% 
  #filter(Solv=="PCM", DFT=="B97D", EnType=="E") %>% 
  rename(x=SD,y=SDLH) %>% 
  mutate(Method=paste0(DFT,"_",Solv,"_",EnType)) %>% 
  group_by(Method) %>%
  do(PredSD=new_r(.,type="continuous")(10**6)) %>% 
  unnest()

write_delim(SDDistrTab,"SDdistributions.csv",";")

SDDistrTab2 <- SDDistrTab %>% 
  group_by(Method) %>%
  mutate(try=row_number())

SDDistrTab3 <- full_join(SDDistrTab2,
                         SDDistrTab2 %>% 
                           rename(Method2=Method,PredSD2=PredSD)) %>% 
  filter(Method2!=Method) %>% 
  mutate(PredSDdiff=PredSD-PredSD2) %>% 
  group_by(Method,Method2) %>% 
  summarize(PredSDmed=median(PredSDdiff),
            PredSD2.5=quantile(PredSDdiff,probs=0.025),
            PredSD97.5=quantile(PredSDdiff,probs=0.975),
            Star=ifelse(PredSD2.5*PredSD97.5<0,"","*"))

SDDistrTab3 %>% filter(Star=="") %>% View()

######### Figure 2D -----

expr <-  data.frame(exp_ee  = c(74, 72, 88, 52, 81, 65, -48, 46, -53, 44,
                                80, 61, -49, 84),
                    temp = c(195, 195, 195, 296, 193, 233, 233, 233, 195,
                             195, 195, 195, 195, 228),
                    Catalyst = c( "1a",   "1b",   "2", "3a",   "3b",   "4a",   "4b",   "4c",   "4d",   "4e",  "5a",   "5b",   "5c",   "6")) 

TestSet <- c("5a","5b","5c","6")

TestTab <- read_tsv("wheeler.csv") %>% 
  melt(id.vars = c("Catalyst", "Product", "TS")) %>%
  rename(Method = variable) %>%
  rowwise() %>% 
  mutate(TS=ifelse(str_detect(Catalyst,"-"),paste0(TS,"_",str_split(Catalyst,"-",simplify=T)[2]),TS),
         Catalyst=str_split(Catalyst,"-",simplify=T)[1]) %>% 
  ungroup() %>% 
  filter(Catalyst%in%TestSet,
         Method=="PCM-B97D_E") %>% 
  group_by(Catalyst,Method) %>% 
  mutate(E=627.5095*(value-min(value,na.rm=T))) %>% 
  inner_join(expr)

TestNaive <- TestTab %>% 
  group_by(Catalyst,Method,exp_ee) %>% 
  summarise(calc_ee=NRJ2ee2(E,Product,temp))

TestNaive %>% 
  mutate(ee_diff=abs(exp_ee-calc_ee))

SDdistr.B97D.PCM.TrSet1 <- AllData %>% 
  filter(Cat=="1a-4e") %>% 
  rename(x=SD,y=SDLH) %>% new_r(type="continuous")

#aaaa <- tibble(a=SDdistr.B97D.PCM.TrSet1(10000))

TestTab_FP_Bay <- TestTab %>% 
  ungroup() %>% 
  group_by(Catalyst, Method) %>% 
  summarise(ee = ForthPropBayes(Product,E,SDdistr.B97D.PCM.TrSet1,temp,100000,F,Catalyst))

B97D.PCM.TrSet1.TestPlotCompar <- TestTab_FP_Bay %>% 
  ggplot(aes(x=ee)) +
  #geom_density(fill=NA,trim = T,alpha=0.7) +
  geom_density_ridges_gradient(aes(y=0,fill=stat(quantile)),quantile_lines = TRUE, 
                               #quantile_fun = hdi, 
                               quantiles =c(0.1,0.9),
                               vline_linetype = 1, key_glyph=draw_key_abline) +
  scale_fill_manual(values = c("transparent", alpha("lightgreen",0.1), "transparent"), guide = "none") +
  #geom_boxplot(fill=NA, ) +
  #scale_color_hue(breaks=c("Bayes","MaxL"),labels=c("Bayesian","Max. Likelihood")) +
  geom_vline(data=TestNaive, 
             aes(xintercept=exp_ee), color="darkgoldenrod1", size=1) +
  geom_vline(data=TestNaive, 
             aes(xintercept=calc_ee), color="#94A5C2", size=1) +
  scale_x_continuous(name=expression(paste("Enantiomeric excess of ",italic("(R)"),", %"))) +
  scale_y_continuous(name="Probability") +
  theme_bw() +
  facet_grid(Catalyst~.,scales="free_y") +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.length.y = unit(0,"mm"),
        strip.text.y = element_text(face = "bold"),
        legend.key.size = unit(2,"mm"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-4,0,0,0,"mm"),
        axis.title.x=element_text(hjust=0.5),
        panel.grid = element_blank())
B97D.PCM.TrSet1.TestPlotCompar
ggsave("B97D.PCM.TrSet1.TestPlotCompar.png",B97D.PCM.TrSet1.TestPlotCompar,width=8.3,height=12.5,dpi=2000, units="cm")

TestTab_FP_Bay %>% group_by(Catalyst) %>% 
  summarise(quant10=quantile(ee,probs = 0.10),
            quant90=quantile(ee,probs = 0.90)) %>% 
  left_join(expr) %>% 
  mutate(right_diff=quant10-exp_ee,
         left_diff=quant90-exp_ee)


##### AUC Plot

AUCtab <- TestTab_FP_Bay %>% 
  mutate(eelev = round(ee,0)) %>% 
  group_by(Catalyst, Method, eelev) %>% 
  summarise(count=n()) %>% 
  arrange(desc(eelev)) %>% 
  mutate(AUC=cumsum(count)/100000) 

AUCtab %>% View()

Lev75pts <- AUCtab %>% filter(eelev==75) %>% 
  select(Catalyst, Method, eelev, AUC)

AUCPlot <- AUCtab %>% 
  ggplot(aes(x=eelev,y=AUC*100, color=Catalyst)) +
  geom_vline(xintercept = 75, color="grey80", size=1.2) +
  geom_line() +
  geom_point(data=Lev75pts) +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  scale_x_continuous(name=expression(paste("Enantiomeric excess of ",italic("(R)")," threshold, %"))) +
  scale_y_continuous(name=expression(paste("Probability of experimental ",italic("ee"),"\nof ", 
                                           italic("(R)"), " being no less, %"))) +
  theme(legend.position = "none",
        #axis.title.x=element_text(hjust=0.5),
        axis.title.y = element_blank(),
        panel.grid = element_blank())
AUCPlot
ggsave("B97D.PCM.TrSet1.AUCPlot2.png",AUCPlot,width=8,height=4,dpi=2000, units="cm")
