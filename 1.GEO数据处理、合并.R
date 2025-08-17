

###加载R包
library(tidyverse)
library(GEOquery)
library(tinyarray)
library(AnnoProbe)
library(limma)
library(stringr)
library(writexl)# 导出数据为xlsx格式
library(data.table)


############################1.GEO数据下载##########################


############################2.读取临床信息##########################
####1.读取数据：GSE135917
###1.读取表达矩阵：
setwd("C:/Users/23634/Desktop/论文代码")
geoID = "GSE135917"
gset = getGEO(geoID, destdir=".", 
              AnnotGPL = F, getGPL = F)#后面的是注释文件，为F的就不下载了


###2.通过pData函数获取分组信息（获取临床数据）
#临床信息数据：
phe=pData(gset[[1]])
table(phe$title) 
#分为两组：
group_list <- ifelse(str_detect(phe$title, "BD"), "BD",
                     "control")#str_detect:字符串中的每个元素都重叠的话返回一个逻辑向量TRUE，这样就不需要把disease state:或者Disease State:列出
table(group_list)
#因子型
group_list = factor(group_list,
                    levels = c("control","BD"))  
table(group_list)
phe$group = group_list#增加一列分组信息
#保存临床信息
write.table(phe,file=paste0(geoID,".clinic.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write.csv(phe, file=paste0(geoID,".clinic.csv"),quote=F)



############################3.提取表达矩阵##########################
#提取表达矩阵
exp <- exprs(gset[[1]])
dim(exp)#看一下dat这个矩阵的维度
exp <- as.data.frame(exp)


############################4.探针ID转换##########################
index = gset[[1]]@annotation#查看平台注释
index
ids <- AnnoProbe::idmap(index)

#对ids的列名进行重命名
colnames(ids) <- c("probe_id", "symbol") #重命名
write.table(ids,file=paste0(geoID,".ids.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write_xlsx(ids, paste0(geoID,".ids.xlsx"))



############################5.基因ID注释##########################
#合并探针和基因名
exp = as.data.frame(exp)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp,by=c("probe_id"="probe_id"))

#去除重复基因
exp2 <- exp1[,-1]
exp <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))



############################6.数据标准化##########################
###以下为limma包对于微矩阵数据的标准化
####看数据是否经过归一化处理，若处理过就不用矫正,直接进行下一步id转换
pdf(file=paste0(geoID,".boxplot-数据标准化前.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=F,col=group_list, las=2)#画图
#标准化之前的counts数据——数据量大且分散
dev.off()
range(exp)#查看数据值范围
write.table(exp,file=paste0(geoID,".数据标准化前.txt"),sep="\t",quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化前.csv"),quote=F)


##（避免多次进行normalizeBetweenArrays）
#exp <- log2(exp+1)#矫正数据
exp=normalizeBetweenArrays(exp)#校正：数据均值已经接近，不需要再次标准化
pdf(file=paste0(geoID,".boxplot-数据标准化后.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)#画图
#标准化之后的vst数据——数据量集中
dev.off()
range(exp)#查看数据值范围

#保存表达量矩阵及分组
save(exp,phe,group_list,file = paste0(geoID,"分组及表达矩阵.Rdata"))
write.table(exp,file=paste0(geoID,".数据标准化后.txt"),sep="\t",,quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化后.csv"),quote=F)

#数据备份：
exp_GSE17114 = exp
phe_GSE17114 = phe
