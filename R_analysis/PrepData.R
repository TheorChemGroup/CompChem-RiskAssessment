library(tidyverse)

alldata <- bind_rows(read_delim("new_data\\PCM-B97D_E_Cats.csv",",") %>% 
  select(-...1, -Method) %>% mutate(Method="PCM.B97D.E"),
  read_delim("new_data\\Distr16.PCM-B97D_E.csv",",") %>% mutate(Cat="All", Method="PCM.B97D.E"),
  read_delim("new_data\\Distr1a-4e.PCM-B97D_E.csv",",") %>% mutate(Cat="1a-4e", Method="PCM.B97D.E"),
  read_delim("new_data\\Distr16.B97D.csv",",") %>% mutate(Cat="All", Method="None.B97D.E"),
  read_delim("new_data\\Distr16.PCM-M06-2X.csv",",") %>% mutate(Cat="All", Method="PCM.M06-2X.E"),
  read_delim("new_data\\Distr16.PCM-wB97XD.csv",",") %>% mutate(Cat="All", Method="PCM.wB97XD.E"),
  read_delim("new_data\\Distr16.PCM-B97D_G.csv",",") %>% mutate(Cat="All", Method="PCM.B97D.G"),
  read_delim("new_data\\Distr16.PCM-B97D_H.csv",",") %>% mutate(Cat="All", Method="PCM.B97D.H"),
  read_delim("new_data\\Distr16.PCM-B97D_qG.csv",",") %>% mutate(Cat="All", Method="PCM.B97D.qhG"),
  read_delim("new_data\\Distr16.SMD-B97D.csv",",") %>% mutate(Cat="All", Method="SMD.B97D.E"),
  read_delim("new_data\\Distr16.SMD-M06-2X.csv",",") %>% mutate(Cat="All", Method="SMD.M06-2X.E"),
  read_delim("new_data\\Distr_random.B97D.csv",",") %>% mutate(Cat="AllRuined", Method="None.B97D.E"),
  read_delim("new_data\\Distr_random.PCM-B97D_E.csv",",") %>% mutate(Cat="AllRuined", Method="PCM.B97D.E"),
  read_delim("new_data\\Distr_random.PCM-B97D_G.csv",",") %>% mutate(Cat="AllRuined", Method="PCM.B97D.G"),
  read_delim("new_data\\Distr_random.PCM-B97D_H.csv",",") %>% mutate(Cat="AllRuined", Method="PCM.B97D.H"),
  read_delim("new_data\\Distr_random.PCM-B97D_qG.csv",",") %>% mutate(Cat="AllRuined", Method="PCM.B97D.qhG"),
  read_delim("new_data\\Distr_random.PCM-M06-2X.csv",",") %>% mutate(Cat="AllRuined", Method="PCM.M06-2X.E"),
  read_delim("new_data\\Distr_random.PCM-wB97XD.csv",",") %>% mutate(Cat="AllRuined", Method="PCM.wB97XD.E"),
  read_delim("new_data\\Distr_random.SMD-B97D.csv",",") %>% mutate(Cat="AllRuined", Method="SMD.B97D.E"),
  read_delim("new_data\\Distr_random.SMD-M06-2X.csv",",") %>% mutate(Cat="AllRuined", Method="SMD.M06-2X.E"))

alldata %>% distinct(Method)

alldata %>% filter(Cat=="1a") %>%
  ggplot(aes(x=Distr)) +
  geom_density(aes(y=after_stat(scaled)), bounds=c(0,Inf))

densfun <- function(distr, npts=24000, adj=1) {
  d = density(distr, kernel="gaussian", from=-1, to=11, n=npts, adjust=adj)
  ddd <- tibble(x=round(d$x,3),y=d$y) %>% 
    mutate(x=ifelse(x>10,20-x,abs(x))) %>% 
    group_by(x) %>% 
    summarize(y=sum(y)) %>%
    ungroup() %>% 
    mutate(y=y/max(y)) %>% 
    filter(x!=0, x!=10)
  return(ddd)
}

alldata %>% filter(Cat=="1a") %>% 
  group_by(Cat,Method) %>% 
  reframe(densfun(Distr)) %>% #View()
  ggplot(aes(x=x,y=y)) +
  geom_line()

alldens <- alldata %>% 
  group_by(Cat,Method) %>% 
  reframe(densfun(Distr)) %>% 
  rename(SD=x, SDLH=y) %>% 
  separate(col=Method,sep="\\.",into=c("Solv","DFT","EnType"))

alldens

alldens %>% write_delim("AllData.csv",";")
