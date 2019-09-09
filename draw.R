require(ggplot2); require(scales); require(reshape2)

#d=(read.csv("allres.csv",sep="\t",header=F))
#d=(read.csv("~/Linux Files/AlignmentErrorRemoval-TestingPipeline/allres.csv", sep="\t", header=F))
#names(d) <- c("E","DR","V3","DR2","V5","Diameter","PD","N","varname","var","V11", "Rep","FP","FN","TP","TN")
d=(read.csv("~/Linux Files/AlignmentErrorRemoval-TestingPipeline/allres_v2.csv", sep=",", header=F))
names(d) <- c("E","DR","V3","DR2","V5","Diameter","PD","N","varname","var","Rep","FP0", "FN0", "TP0", "TN0", "FP","FN","TP","TN")

nlabels = c("1","2%","5%","10%","20%")

# 16S.B: K - Recall vs Diameter
ggplot(aes(x=Diameter,y=TP/(TP+FN), 
      color=as.factor(var)),data=d[d$E=="16S.B_K",])+geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="k")
ggsave("Recall-k.pdf",width = 6,height = 6)


# 16S.B: K - FPR vs Diameter
ggplot(aes(x=Diameter,y=FP/(FP+TN)+0.0001, 
      color=as.factor(var)),data=d[d$E=="16S.B_K" & d$var<39,])+geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="k")+scale_y_log10(name="FPR")
ggsave("FPR-k.pdf",width = 6,height = 6)


# 16S.B: K - Recall vs K (violin)
ggplot(aes(x=var,y=TP/(TP+FN), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1)),
      color=Diameter),data=d[d$E=="16S.B_K" & d$var<45,])+geom_jitter(alpha=0.4)+geom_violin(aes(group=var),fill=NA)+
  theme_classic()+scale_x_continuous(breaks = c(5,11,17,23,29,35,41),name="k")+
  scale_y_continuous(name="Recall")+scale_shape_manual(name="FPR",values=c(1,16,4))
ggsave("violion-k.pdf",width = 6,height = 6)


# 16S.B: K - Recall vs FPR
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var)),data=d[d$E=="16S.B_K" & d$var<42,])+
  geom_point(alpha=0.25)+geom_smooth()+
  theme_classic()+theme(legend.position = c(.85,.25))+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="k")+coord_cartesian(xlim=c(0,0.003))+
  scale_x_continuous(name="FPR")+scale_y_continuous("Recall")
ggsave("ROC-k.pdf",width = 6,height = 6)


# 16S.B: K - FPR vs Recall
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      shape=as.factor(var), linetype=as.factor(var),colour=Diameter),data=d[d$E=="16S.B_K" & d$var<39,])+
  geom_point(alpha=0.5)+geom_smooth(se=F,color="black")+
  theme_classic()+#theme(legend.position = c(.85,.25))+
  scale_shape(name="k")+scale_linetype(name="k")+
  coord_cartesian(xlim=c(0,0.003))+
  scale_x_continuous(name="Recall")+scale_y_continuous("FPR")
ggsave("ROC-k-2.pdf",width = 6,height = 6)


# 16S.B: ErrLen - Recall vs Error Length (Violin)
ggplot(aes(x=var,y=TP/(TP+FN), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1)),
      color=Diameter),data=d[d$E=="16S.B_ErrLen" ,])+geom_jitter(alpha=0.3)+geom_violin(aes(group=var),fill=NA)+
  theme_classic()+scale_x_log10(breaks=c(2,3,4,8,16,32,64),name="error length (x times k)")+
  scale_y_continuous(name="Recall")+scale_shape_manual(name="FPR",values=c(1,16,4,3))
ggsave("violion-len.pdf",width = 6,height = 6)


# 16S.B: ErrLen - FPR vs Diameter
ggplot(aes(x=Diameter,y=FP/(FP+TN)+0.0001,color=as.factor(var)),data=d[d$E=="16S.B_ErrLen",])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len")+
  scale_y_log10(name="FPR")+
  geom_vline(aes(xintercept = (FP0!=0)*Diameter), alpha=0.1)

ggsave("FPR-errlen.pdf",width = 6,height = 6)


# 16S.B: ErrLen - Recall vs Diameter 
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(var)),data=d[d$E=="16S.B_ErrLen",])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "k (=11)")))+
  geom_vline(aes(xintercept = (FP0!=0)*Diameter), alpha=0.1)
ggsave("Recall-errlen.pdf",width = 6,height = 6)


# 16S.B: NumErrSeq - Recall vs Number of Sequences with Error (violin)
ggplot(aes(x=as.factor(100/((var!=N)*var+(var==N)*100)),group=as.factor(100/((var!=N)*var+(var==N)*100)),y=TP/(TP+FN), 
      color=cut(Diameter,breaks=5), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1))),
      data=d[d$E=="16S.B_NumErrSeq" ,])+
  geom_jitter(alpha=0.4)+geom_violin(fill=NA)+
  theme_classic()+scale_x_discrete(name="Number of sequences with error",labels=nlabels)+
  scale_y_continuous(name="Recall")+scale_shape_manual(name="FPR",values=c(1,16,4))+
  scale_color_brewer(name="Diameter",palette = "Dark2")
ggsave("violion-N.pdf",width = 6,height = 6)


# 16S.B: NumErrSeq: FPR vs Diameter
ggplot(aes(x=Diameter,y=FP/(FP+TN)+0.00001,color=as.factor(100/((var!=N)*var+(var==N)*100))),data=d[d$E=="16S.B_ErrSeq" ,])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="n")+
  scale_y_log10(name="FPR")
ggsave("FPR-N.pdf",width = 6,height = 6)


# 16S.B: NumErrSeq: Recall vs Diameter
ggplot(aes(x=Diameter,y=TP/(TP+FN), group= as.factor(100/((var!=N)*var+(var==N)*100)),
      color=as.factor(100/((var!=N)*var+(var==N)*100)), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1))),
      data=d[d$E=="16S.B_ErrSeq" ,])+geom_point(alpha=0.5)+
  theme_classic()+geom_smooth(se=F)+scale_shape_manual(name="FPR",values=c(1,16,4))+
  scale_color_brewer(palette = "Paired",name="n",labels=nlabels)+
  scale_y_continuous(name="Recall")+coord_cartesian(ylim=c(0.35,1))+
  geom_vline(aes(xintercept = (FP0!=0)*Diameter), alpha=0.1)
ggsave("Recall-N.pdf",width = 6,height = 6)


# 16S.B: NumErrSeq: Recall vs FPR
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN),color=as.factor(100/((var!=N)*var+(var==N)*100))),data=d[d$E=="16S.B_NumErrSeq",])+
  geom_point(alpha=0.25)+#geom_smooth(se=F)+
  theme_classic()+theme(legend.position = c(.85,.25))+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="n",labels=nlabels)+coord_cartesian(xlim=c(0,0.003))+
  scale_x_continuous(name="FPR")+scale_y_continuous("Recall")
ggsave("ROC-n.pdf",width = 6,height = 6)


# 16S.B: NumErrSeq: FPR vs Recall
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN),shape=as.factor(100/((var!=N)*var+(var==N)*100)),colour=Diameter),data=d[d$E=="16S.B_NumErrSeq",])+
  geom_point(alpha=0.5)+
  theme_classic()+#theme(legend.position = c(.85,.25))+
  scale_shape(name="n",labels=nlabels)+scale_linetype(name="n",labels=nlabels)+
  coord_cartesian(xlim=c(0,0.003))+
  scale_x_continuous(name="Recall")+scale_y_continuous("FPR")
ggsave("ROC-n-2.pdf",width = 6,height = 6)


# Hackett: ErrLen - Recall vs Diameter 
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(var)),data=d[d$E=="Hackett_ErrLen" & d$var<32,])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "k (=11)")))+
  geom_vline(aes(xintercept = (FP0!=0)*Diameter), alpha=0.1)
ggsave("Hackett-RecallvDiameter-ErrLen.pdf",width=6,height=6)


# Hackett: K - Recall vs Diameter 
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(var)),data=d[d$E=="Hackett_K",])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="k")
ggsave("Hackett-RecallvDiameter-K.pdf",width=6,height=6)


# AA-RV100-BBA0039: ErrLen - Recall vs Diameter 
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(var)),data=d[d$E=="small-10-aa-RV100-BBA0039_ErrLen" & d$var<32,])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "k (=11)")))+
  geom_vline(aes(xintercept = (FP0!=0)*Diameter), alpha=0.1)
ggsave("AA-RecallvDiameter-ErrLen.pdf",width=6,height=6)
  

# AA-RV100-BBA0039: K - Recall vs Diameter 
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(var)),data=d[d$E=="small-10-aa-RV100-BBA0039_K",])+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="k")
ggsave("AA-RecallvDiameter-K.pdf",width=6,height=6)





summ_roc <- function(d2,form) {
  ad2 = dcast(d2, form ,fun.aggregate=sum,value.var = c("FP"))
  ad2=cbind(dcast(d2, form ,fun.aggregate=sum,value.var = c("FP")),
            dcast(d2, form ,fun.aggregate=sum,value.var = c("FN"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("TP"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("TN"))[,length(ad2)]
            )
  names(ad2)[(length(names(ad2))-3):(length(names(ad2)))]=c("FP","FN","TP","TN")
  ad2
}


# 16S.B: NumErrSeq - Recall vs FPR (sum)
d2=d[d$E=="16S.B_ErrSeq",]
d2$n=with(d2,as.factor(100/((var!=N)*var+(var==N)*100)))
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), color=as.factor(n) ),data=summ_roc(d2,n~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = c(.85,.25))+
  scale_color_brewer(name="n",labels=nlabels, palette="Paired")+
  scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")
ggsave("sum-N.pdf",width=4,height = 4)


# 16S.B: K - Recall vs FPR (sum)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var)),data=summ_roc(d[d$E=="16S.B_K" & d$var<45,],var~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = c(.85,.35))+
  scale_shape(name="k")+scale_color_brewer(name="k",palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")
ggsave("sum-k.pdf",width=4,height = 4)


# 16S.B: K - Recall vs FPR (sum) [includes k and diameter] 
options(digits = 2)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
           color=as.factor(var),shape=`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`),data=summ_roc(d[d$E=="16S.B_K",],                                                                                                         var+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Diameter")+scale_color_brewer(name="k",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.13,1))
ggsave("Hackett-ROC-K.pdf",width=6,height = 6)


# 16S.B: ErrLen - Recall vs FPR (sum)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var)),data=summ_roc(d[d$E=="16S.B_ErrLen",],var~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = c(.85,.35))+
  scale_shape(name="k")+scale_color_brewer(name="error len",palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")
ggsave("sum-len.pdf",width=4,height = 4)


# 16S.B: ErrLen - Recall vs FPR (sum) [includes error length and diameter] 
options(digits = 2)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var),shape=`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`),data=summ_roc(d[d$E=="16S.B_ErrLen",],
      var+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error len",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "k (=11)")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.13,1))
ggsave("sum-len-diam.pdf",width=5,height = 4.2)


# 16S.B: ErrLen - Recall vs FPR (sum) [includes error length and percentage of sequences with error]
ad3 =  summ_roc(d2,n~.) 
ad3$var=8
ad3 = rbind(ad3,data.frame(n=5,summ_roc(d[d$E=="16S.B_ErrLen",],var~.))[,c(1,3:6,2)])
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN),color=as.factor(var),shape=as.factor(n)),data=ad3)+
  geom_point(alpha=1)+
  theme_light()+theme()+
  scale_shape(name="n",labels=nlabels)+
  scale_color_brewer(name="error len",palette = "Paired",labels = function(x) (paste(x,"k")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.2,1))
ggsave("sum-len-n.pdf",width=4.5,height = 4)


# Hackett: K - Recall vs FPR (sum) [includes k and diameter] 
options(digits = 2)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var),shape=`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`),data=summ_roc(d[d$E=="Hackett_K",],
      var+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Diameter")+scale_color_brewer(name="k",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.13,1))
ggsave("Hackett-ROC-K.pdf",width=6,height = 6)


# Hackett: ErrLen - Recall vs FPR (sum) [includes k and diameter] 
options(digits = 2)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var),shape=`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`),data=summ_roc(d[d$E=="Hackett_ErrLen",],
      var+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "k (=11)")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.13,1))
ggsave("Hackett-ROC-ErrLen.pdf",width=6,height=6)


# AA-RV100-BBA0039: K - Recall vs FPR (sum) [includes k and diameter] 
options(digits = 2)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var),shape=`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`),data=summ_roc(d[d$E=="small-10-aa-RV100-BBA0039_K",],
      var+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Diameter")+scale_color_brewer(name="k",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.13,1))
ggsave("AA-ROC-K.pdf",width=6,height=6)


# AA-RV100-BBA0039: ErrLen - Recall vs FPR (sum) [includes k and diameter] 
options(digits = 2)
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), 
      color=as.factor(var),shape=`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`),data=summ_roc(d[d$E=="small-10-aa-RV100-BBA0039_ErrLen",],
      var+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "k (=11)")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+coord_cartesian(ylim=c(0.13,1))
ggsave("AA-ROC-ErrLen.pdf",width=6,height=6)
