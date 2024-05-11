library(readr);library(magrittr);library(dplyr);library(DescTools);library(MASS);library(randomForestSRC)
library(ggplot2);library(survival);library(caret);library(keras);library(survminer);library(readxl)
library(timeROC)
`%>%` <- magrittr::`%>%` 
data <- read_excel("SEER3.xlsx")
#dd <- data[data$`Year of diagnosis`==2016|data$`Year of diagnosis`==2017,]
#write.csv(dd,file = "1617.csv")
#data <- data[data$`Year of diagnosis` >=2000 & data$`Year of diagnosis`<= 2015,]
data <- data[,c(1,2,6,3:5,16:18,24,28,29,37,41,42,51,58:63,45,46,50,55,56,57,64,65)]
##############stage Ⅲ##########################
#####EOD###############
# 检查数据框的列名是否有重复
duplicate_names <- names(data)[duplicated(names(data))]
if(length(duplicate_names) > 0) {
  cat("发现重复列名: ", duplicate_names, "\n")
  # 使列名唯一
  names(data) <- make.unique(names(data))
} else {
  cat("没有重复列名。\n")
}
#####2004-2015
#data1 <- data %>%
#  mutate(stage = ifelse(`CS extension (2004-2015)` %in% c(130, 140, 150, 160, 170, 200, 300, 400, 410, 420, 450, 458, 460, 470, 500, 550, 560, 565, 570, 600, 650, 655, 660, 675, 700, 750, 800, 850, 900) &
#                          `CS lymph nodes (2004-2015)` %in% c(050, 100, 110, 200, 210, 220, 300, 400, 410, 420, 430, 450, 460, 470, 480, 800) &
 #                         `CS mets at dx (2004-2015)` == "0","3期", "0"))

#data1 <- data1[data1$stage=="3期",]
########2004-2015
data1 <- data %>%
 filter(`Derived AJCC Stage Group, 6th ed (2004-2015)` %in% c("IIIA","IIIB", "IIIC", "IIINOS")
        )
#######2016-2017
data2 <- data %>%
  filter(`7th Edition Stage Group Recode (2016-2017)` %in% c("III","IIIA", "IIIB", "IIIC"))
######2018+
data3 <- data %>%
  filter(`Derived EOD 2018 Stage Group (2018+)` %in% c("3","3A", "3B", "3C"))
data <- rbind(data1,data2,data3)
data <- data[data$`Survival months` !="Unknown",]   
data$`Survival months` <- as.numeric(data$`Survival months`)
data <- data[data$`Race recode (White, Black, Other)`!= "Unknown",]
d <- data[data$`Survival months` != 0,]
d <- d[!(d$`Vital status recode (study cutoff used)`=="Alive"& d$`Survival months`<60),]

#年龄分组
#table(d$`Age recode with <1 year olds`)
# 提取年龄
d$age1 <- as.numeric(sub("([0-9]+)(-.+|\\+.+)", "\\1", d$`Age recode with <1 year olds`))

# 使用 cut 函数创建年龄分组
d$agegroup <- cut(d$age1, 
                  breaks = c(-Inf, 44, 69, Inf), 
                  labels = c("1", "2", "3"),
                  right = FALSE)
d <- d[!is.na(d$agegroup), ]
#################婚姻状态分组#############################
d <- d[d$`Marital status at diagnosis`!="Unknown",]
d$`Marital status at diagnosis`[d$`Marital status at diagnosis`=="Married (including common law)"] <- "Married"
d$`Marital status at diagnosis`[d$`Marital status at diagnosis`=="Single (never married)"] <- "Single"
d$`Marital status at diagnosis`[d$`Marital status at diagnosis` !="Single" & d$`Marital status at diagnosis` !="Married"] <- "SDM"
table(d$`Marital status at diagnosis`)
sum(is.na(d$`Marital status at diagnosis`))

table(d$`Median household income inflation adj to 2021`)
##################income###############################
df <- d %>%
  mutate(
    income = case_when(
      d$`Median household income inflation adj to 2021`  %in% c("< $35,000","$35,000 - $39,999", "$40,000 - $44,999", "$45,000 - $49,999", "$50,000 - $54,999") ~ "Group1",
      d$`Median household income inflation adj to 2021`  %in% c("$55,000 - $59,999", "$60,000 - $64,999", "$65,000 - $69,999", "$70,000 - $74,999") ~ "Group2",
      d$`Median household income inflation adj to 2021`  %in% c("$75,000+") ~ "Group3",
      TRUE ~ NA_character_  # 如果没有匹配任何组，可以分配 NA 或者其他默认值
    )
  )
df$income
table(df$income)
sum(is.na(df$income))

#########################rural urban####################
df <- df %>%
  mutate(
    Continuum = case_when(
      df$`Rural-Urban Continuum Code`  %in% c("Counties in metropolitan areas ge 1 million pop", "Counties in metropolitan areas of 250,000 to 1 million pop", "Counties in metropolitan areas of lt 250 thousand pop") ~ "Counties",
      df$`Rural-Urban Continuum Code`  %in% c("Nonmetropolitan counties adjacent to a metropolitan area", "Nonmetropolitan counties not adjacent to a metropolitan area") ~ "Nonmetropolitan counties",
      TRUE ~ NA_character_  # 如果没有匹配任何组，可以分配 NA 或者其他默认值
    )
  )
table(df$Continuum)
sum(is.na(df$Continuum))
df <- df[!is.na(df$Continuum),]

#################radiation############
table(df$`Radiation recode`)
df$radiation <- ifelse(df$`Radiation recode` == "None/Unknown", 0, 1)
sum(is.na(df$radiation))

###############chemo##################
table(df$`Chemotherapy recode (yes, no/unk)`)
df$chemo <- df$`Chemotherapy recode (yes, no/unk)`
sum(is.na(df$chemo))

#########################surgery seq#########
table(df$`RX Summ--Surg/Rad Seq`)
df$rad_surg_sequence <- with(df, ifelse(`RX Summ--Surg/Rad Seq` == "No radiation and/or cancer-directed surgery", "No Radiation",
                                        ifelse(`RX Summ--Surg/Rad Seq` == "Radiation prior to surgery", "Radiation Before Surgery",
                                               ifelse(`RX Summ--Surg/Rad Seq` == "Radiation after surgery", "Radiation After Surgery",
                                                      ifelse(`RX Summ--Surg/Rad Seq` == "Radiation before and after surgery", "Radiation Before and After Surgery",
                                                             "Other")))))

df$rad_surg_sequence <- as.factor(df$rad_surg_sequence)
#df <- na.omit(df)
df <- df[!df$`Total number of in situ/malignant tumors for patient`=="Unknown",]
sum(is.na(df$`Total number of in situ/malignant tumors for patient`))
df <- df[!is.na(df$`Total number of in situ/malignant tumors for patient`),]

########基线特征#############
count_table <- table(df$agegroup)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

count_table <- table(df$Sex)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

count_table <- table(df$`Race recode (White, Black, Other)`)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

count_table <- table(df$`Marital status at diagnosis`)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

df$numberT <- as.numeric(df$`Total number of in situ/malignant tumors for patient`)
df$numberT[df$numberT>=3] <- "3"
count_table <- table(df$numberT)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table
df$`Total number of in situ/malignant tumors for patient` <- as.numeric(df$`Total number of in situ/malignant tumors for patient`)
summary(df$`Total number of in situ/malignant tumors for patient`)

count_table <- table(df$`Chemotherapy recode (yes, no/unk)`)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

count_table <- table(df$income)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

count_table <- table(df$Continuum)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

df$status <- df$`Vital status recode (study cutoff used)`
df$status[df$`SEER other cause of death classification`=="Alive or dead due to cancer" & df$`Vital status recode (study cutoff used)`=="Dead"] <- "2"
df$status[df$status != "2"] <- "1"
table(df$status)
count_table <- table(df$status)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

df$histology <- df$`Histology recode - broad groupings`
df$histology[df$histology=="8140-8389: adenomas and adenocarcinomas"] <- "Adenocarcinomas"
df$histology[df$histology=="8440-8499: cystic, mucinous and serous neoplasms"] <- "CMS neoplasms"
df$histology[df$histology !="CMS neoplasms" & df$histology !="Adenocarcinomas"] <- "Others"
table(df$histology)
count_table <- table(df$histology)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

table(df$`Regional nodes examined (1988+)`)
df$`Regional nodes examined (1988+)` <- as.numeric(df$`Regional nodes examined (1988+)`)
summary(df$`Regional nodes examined (1988+)`)

table(df$radiation)
count_table <- table(df$radiation)
percentage_table <- prop.table(count_table) * 100
count_table
percentage_table

table(df$rad_surg_sequence)

table(df$`Regional nodes positive (1988+)`)
df$`Regional nodes positive (1988+)` <- as.numeric(df$`Regional nodes positive (1988+)`)
summary(df$`Regional nodes positive (1988+)`)

table(df$`Survival months`)
summary(df$`Survival months`)

table(df$`RX Summ--Surg/Rad Seq`)


#############################差异性检验#######################
x<-xtabs(~df$Sex+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)     

x<-xtabs(~df$agegroup+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)    

x<-xtabs(~df$`Marital status at diagnosis`+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)    

x<-xtabs(~df$`Race recode (White, Black, Other)`+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)    

x<-xtabs(~df$numberT+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)   

x<-xtabs(~df$`Chemotherapy recode (yes, no/unk)`+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F) 

x<-xtabs(~df$radiation+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)

x<-xtabs(~df$histology+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F) 

x<-xtabs(~df$radiation+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)

x<-xtabs(~df$chemo+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)

x<-xtabs(~df$income+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)

x<-xtabs(~df$rad_surg_sequence+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)

x<-xtabs(~df$Continuum+df$status,data=df)
addmargins(x)
addmargins(prop.table(x, 1)*100, 2)
chisq.test(x,correct=F)

summary(df$`Regional nodes examined (1988+)`)
wilcox.test(df$`Regional nodes examined (1988+)`~df$status,data=df,exact=F) #（有统计学意义）
wilcox.test(df$`Regional nodes positive (1988+)` ~df$status,data=df,exact=F) #（有统计学意义）
wilcox.test(df$`Survival months` ~df$status,data=df,exact=F) #（有统计学意义）
wilcox.test(df$`Total number of in situ/malignant tumors for patient` ~df$status,data=df,exact=F) #（有统计学意义）


####################cochran趋势性检验###########################
# 仅选择化疗方案记录为'yes'的数据
# 筛选出2023年的数据
subset_df = df[df$`Year of diagnosis` == 2017, ]
# 计算其中Chemotherapy recode为"Yes"的数量
yes_count = table(subset_df$`Chemotherapy recode (yes, no/unk)`)["Yes"]


####设置罗马字体#####
#windowsFonts(roman = windowsFont("TT Times New RomanS"))
par(family = "roman")
df_summary <- df %>%
  group_by(`Year of diagnosis`) %>%
  summarise(
    Total = n(),  # 计算每年的总患者数
    YesCount = sum(`Chemotherapy recode (yes, no/unk)` == "Yes", na.rm = TRUE),  # 计算每年选择"yes"的患者数
    YesPercentage = round(YesCount / Total * 100,2)  # 计算百分比
  ) %>%
  ungroup()  # 移除分组
ggplot(data = df_summary, aes(x = `Year of diagnosis`, y = YesPercentage)) +
  geom_line(aes(color = "Line"), size = 1) +  # 将颜色属性添加到aes内，并给予名字"Line"
  geom_point(aes(color = "Point"), size = 3) + # 同样将颜色属性添加到aes内，并命名为"Point"
  geom_text(aes(label = YesPercentage), vjust = -1, color = "black", size = 3) + 
  scale_color_manual(values = c("Line" = "skyblue"#, "Point" = "black"
                                ), 
                     name = "", 
                     labels = c("Trend Line"#, "Patients with guideline-recommended chemotherapy (100%)"
                                )) + # 自定义图例的颜色和标签
  scale_x_continuous(breaks = 2000:2020) + 
  theme_minimal() +  
  labs(x = "Year of Diagnosis", y = "Percentage of Patients (100%)")+
  ggtitle("The Percentage of Patients Receiving Guideline-Recommended Chemotherapy Over the Years")+
       #title = "The Percentage of Patients Receiving Guideline-Recommended Chemotherapy Over the Years") +
  theme(plot.title = element_text(hjust = 0.2,size = 9, face = "bold"))
df$TreatmentBinary <- ifelse(df$`Chemotherapy recode (yes, no/unk)` == "Yes", 1, 0)
#使用Cochran-Armitage趋势检验来评估趋势
tab <- table(df$`Year of diagnosis`, df$TreatmentBinary)
result <- CochranArmitageTest(tab)
print(result)

###############################overall survival######
survobj <- data.frame(df$`Survival months`, df$status) # 这里主要是指定时间和生存状态
colnames(survobj) <- c("time","status")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ 1, data = survobj)
# 使用 plot() 绘制生存曲线
plot(fit0, xlab = "Survival Time (Months)", ylab = "Survival Probability", main = "Survival Curve")
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv[1] * 100, "%"))
print(paste("5年生存率：", surv_5years$surv[1] * 100, "%"))
surv_3years
surv_5years

##############年龄分组##################
survobj$age <- df$agegroup
colnames(survobj) <- c("time","status","age")
fit0 <- survfit(Surv(time,status) ~ age, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Age)",
  legend.labs = c("<45","[45,69)","≥69"), 
  legend.title = "Age",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ age, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

############性别#################
survobj <- data.frame(df$`Survival months`, df$status,df$Sex)
colnames(survobj) <- c("time","status","sex")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ sex, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Gender)",
  legend.labs = c("Female","Male"), 
  legend.title = "Gender",  
  xlim = c(0, 200),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ sex, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

################Race#################
survobj <- data.frame(df$`Survival months`, df$status,df$`Race recode (White, Black, Other)`)
colnames(survobj) <- c("time","status","race")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ race, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Race)",
  legend.labs = c("Black ","Other","White"), 
  legend.title = "Race",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ race, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)


################Marital#################
survobj <- data.frame(df$`Survival months`, df$status,df$`Marital status at diagnosis`)
colnames(survobj) <- c("time","status","marital")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ marital, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Marital)",
  legend.labs = c("Married ","SDM","Single"), 
  legend.title = "Marital",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标http://127.0.0.1:39067/graphics/9e8508d8-af0c-47af-a38f-8b41b1d2d313.png题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ marital, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)


################Chemotherapy#################
survobj <- data.frame(df$`Survival months`, df$status,df$`Chemotherapy recode (yes, no/unk)`)
colnames(survobj) <- c("time","status","chemo")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ chemo, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Chemotherapy)",
  legend.labs = c("No/Unknown","Yes"), 
  legend.title = "Chemotherapy",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ chemo, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

################Radiation#################
survobj <- data.frame(df$`Survival months`, df$status,df$radiation)
colnames(survobj) <- c("time","status","radiation")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ radiation, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Radiation)",
  legend.labs = c("No/Unknown","Yes"), 
  legend.title = "Radiation",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ radiation, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

################Histology#################
survobj <- data.frame(df$`Survival months`, df$status,df$histology)
colnames(survobj) <- c("time","status","histo")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ histo, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (histology)",
  legend.labs = c("Adenocarcinomas","CMS neoplasms","Others"), 
  legend.title = "Histology",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ histo, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

#################number of tumor###############
survobj <- data.frame(df$`Survival months`, df$status,df$numberT)
colnames(survobj) <- c("time","status","numberT")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ numberT, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (numberT)",
  legend.labs = c("1","2","3+"), 
  legend.title = "numberT",  
  xlim = c(0, 200),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ numberT, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

#################income###############
survobj <- data.frame(df$`Survival months`, df$status,df$income)
colnames(survobj) <- c("time","status","income")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ income, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (income)",
  legend.labs = c("<$35,000-$54,999","$55,000-$74,999",">$74,999"), 
  legend.title = "Income",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ income, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

#################Continuum###############
survobj <- data.frame(df$`Survival months`, df$status,df$Continuum)
colnames(survobj) <- c("time","status","Continuum")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ Continuum, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (Continuum)",
  legend.labs = c("Metropolitan Counties","Nonmetropolitan counties"), 
  legend.title = "Continuum",  
  xlim = c(0, 135),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ Continuum, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

##################rad_chemo_seq#############
survobj <- data.frame(df$`Survival months`, df$status,df$rad_surg_sequence)
colnames(survobj) <- c("time","status","rad_surg_seq")
survobj$status <- as.numeric(survobj$status)
fit0 <- survfit(Surv(time,status) ~ rad_surg_seq, data = survobj)
# 获取3年和5年的生存概率
surv_3and5years <- summary(fit0, times = c(36, 60))
summary(fit0)$table
# 生成生存曲线
plot <- ggsurvplot(
  fit0,                     
  data = survobj,           
  risk.table = TRUE,        
  pval = TRUE,              
  conf.int = TRUE,          
  palette = "jco",          
  xlab = "Survival Time (Months)",        
  ylab = "Survival Probability", 
  title = "Survival Curve (radiation sequence)",
  legend.labs = c("No Ra","Other","Ra-aft-sur","Ra-b/a-sur","Ra-bef-sur"), 
  legend.title = "Radiation sequnce",  
  xlim = c(0, 200),         
  risk.table.height = 0.3  # 调整number at risk表格的高度
)
# 修改生存曲线的ggplot对象
plot$plot <- plot$plot + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2)  # 调整标题位置到图形的正上方
  )
print(plot$plot)
# 获取3年生存率和5年生存率的估计值
surv_3years <- summary(fit0, times = c(36))
surv_5years <- summary(fit0, times = c(60))
# 输出结果
print(paste("3年生存率：", surv_3years$surv * 100, "%"))
print(paste("5年生存率：", surv_5years$surv * 100, "%"))
log_rank_test <- survdiff(Surv(time,status) ~ rad_surg_seq, data = survobj)
# 打印log-rank检验结果
print(log_rank_test)

####################cox####################
df$out <- as.numeric(df$status)-1
surv_object <- Surv(time = df$`Survival months`, event = df$out)
df_clean <- data.frame(
  status = df$out,
  #time = df$`Survival months`,
  Age = df$agegroup,
  Histology = df$histology,
  NumberT = df$`Total number of in situ/malignant tumors for patient`,
  Continuum = df$Continuum,
  Income = df$income,
  MaritalStatus = df$`Marital status at diagnosis`,
  NodesExamined = df$`Regional nodes examined (1988+)`,
  NodesPositive = df$`Regional nodes positive (1988+)`,
  Chemotherapy = df$`Chemotherapy recode (yes, no/unk)`,
  Sex = df$Sex,
  Race = df$`Race recode (White, Black, Other)`,
  Radiation=df$radiation
)
# 创建Cox比例风险模型
model <- coxph(surv_object ~ as.factor(Age)+as.factor(Histology)+NumberT+as.factor(Continuum)+
                 as.factor(Income)+as.factor(MaritalStatus)+NodesExamined+
                 NodesPositive+Chemotherapy+
                 Sex+as.factor(Race)+Radiation, data = df_clean[,-1])

# 查看模型摘要
model_summary <- summary(model)
# 计算优势比 (OR) 和相应的 95% 置信区间 (CI)
or <- exp(coef(model))
ci <- exp(confint(model))
# 将OR和CI整理成一个表格
result_table <- data.frame(
  Variable = names(or),
  OR = or,
  LowerCI95 = ci[,1],
  UpperCI95 = ci[,2],
  P_Value = model_summary$coefficients[, "Pr(>|z|)"]
)
# 查看结果
result_table
write.csv(result_table,"multi cox.csv")

cox.zph(model)
fit <- survfit(model)
plot(fit, xlab = "Survival Time (Months)", ylab = "Survival Probability", main = "Survival Curve")

#################multivariate logistic  status在这里有根据时间进行转换#################
df_clean <- data.frame(
  status = df$out,
  #time = df$`Survival months`,
  Age = df$agegroup,
  Histology = df$histology,
  NumberT = df$`Total number of in situ/malignant tumors for patient`,
  Continuum = df$Continuum,
  Income = df$income,
  MaritalStatus = df$`Marital status at diagnosis`,
  NodesExamined = df$`Regional nodes examined (1988+)`,
  NodesPositive = df$`Regional nodes positive (1988+)`,
  Chemotherapy = df$`Chemotherapy recode (yes, no/unk)`,
  Sex = df$Sex,
  Race = df$`Race recode (White, Black, Other)`,
  Radiation=df$radiation,
  time=df$`Survival months`
)
df_clean$out <- as.factor(with(df, ifelse(df$`Vital status recode (study cutoff used)` == "Alive", 0, 1)))
model <- glm(out~as.factor(Age)+as.factor(Histology)+NumberT+as.factor(Continuum)+
               as.factor(Income)+as.factor(MaritalStatus)+NodesExamined+
               NodesPositive+Chemotherapy+
               Sex+as.factor(Race)+Radiation+time, data = df_clean, family = "binomial")
# 进行逐步选择变量
step_model <- stepAIC(model, direction = "both")
# 获取模型的摘要
model_summary <- summary(step_model)
AIC(step_model)
# 计算优势比 (OR) 和相应的 95% 置信区间 (CI)
or <- exp(coef(step_model))
ci <- exp(confint(step_model))
# 将OR和CI整理成一个表格
result_table <- data.frame(
  Variable = names(or),
  OR = or,
  LowerCI95 = ci[,1],
  UpperCI95 = ci[,2],
  P_Value = model_summary$coefficients[, "Pr(>|z|)"]
)
# 查看结果
result_table
write.csv(result_table,"multi logistic.csv")

##########rf survival#########################################
#建立单纯的随机生存森林模型
library(caret)
library(randomForest)
library(pROC)
# Train the Random Forest model with 5-fold cross-validation
set.seed(123) # for reproducibility
summary(df$out)
# 为所有因子变量的水平生成有效的变量名
clean_factor_levels <- function(factor_column) {
  levels(factor_column) <- make.names(levels(factor_column))
  return(factor_column)
}
# 应用这个函数到所有因子类型的变量
df$agegroup <- clean_factor_levels(df$agegroup)
df$Continuum <- clean_factor_levels(df$Continuum)
df$`Marital status at diagnosis` <- clean_factor_levels(df$`Marital status at diagnosis`)
df$`Chemotherapy recode (yes, no/unk)` <- clean_factor_levels(df$`Chemotherapy recode (yes, no/unk)`)
df$Sex <- clean_factor_levels(df$Sex)
df$`Race recode (White, Black, Other)` <- clean_factor_levels(df$`Race recode (White, Black, Other)`)
df$radiation <- clean_factor_levels(df$radiation)
df$income <- clean_factor_levels(df$income)
df$histology <- clean_factor_levels(df$histology)
df$`Total number of in situ/malignant tumors for patient` <- clean_factor_levels(df$`Total number of in situ/malignant tumors for patient`)
df$`Regional nodes examined (1988+)` <- clean_factor_levels(df$`Regional nodes examined (1988+)`)
df$`Regional nodes positive (1988+)` <- clean_factor_levels(df$`Regional nodes positive (1988+)`)
df$time <- clean_factor_levels(df$`Survival months`)
##########筛选变量 这里并没有筛选变量 全部纳入#############
df_clean <- data.frame(
  status = df$out,
  time = df$`Survival months`,
  Age = df$agegroup,
  Histology = df$histology,
  NumberT = df$`Total number of in situ/malignant tumors for patient`,
  Continuum = df$Continuum,
  Income = df$income,
  MaritalStatus = df$`Marital status at diagnosis`,
  NodesExamined = df$`Regional nodes examined (1988+)`,
  NodesPositive = df$`Regional nodes positive (1988+)`,
  Chemotherapy = df$`Chemotherapy recode (yes, no/unk)`,
  Sex = df$Sex,
  Race = df$`Race recode (White, Black, Other)`,
  Radiation=df$radiation
  #Time=df$time
)
library(randomForest)
#result = replicate(5, rfcv(trainx = df_clean[,-1],
#                           trainy = df_clean$status,
#                           cv.fold = 5),
#                           simplify = FALSE)
#error.cv=sapply(result,"[[","error.cv")
#matplot(result[[1]]$n.var,cbind(rowMeans(error.cv),error.cv),type = "l",lwd = c(2,rep(1,ncol(error.cv))),col = 1,lty = 1,log="x",xlab = "Number of variables",ylab="CV Error")
#提取验证结果绘图
#otu_train.cv <- data.frame(sapply(result, '[[', 'error.cv'))
#otu_train.cv <- sqrt(otu_train.cv)
#otu_train.cv$otus <- rownames(otu_train.cv)
#otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
#otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

#otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
#head(otu_train.cv.mean, 10)
#ggplot(otu_train.cv.mean, aes(Group.1, x)) +
#  geom_line() +
#  theme(axis.title.x =element_text(size=13), axis.title.y=element_text(size=13))+
#  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
#  labs(title = '',x = "Number of variables", y = 'Cross-validation mean RMSE',size=1)
#ggsave("RFvariable.tiff", units="in", width=10, height=9, dpi=500, compression = 'lzw')

###变量重要性
#Example code for Random Forest model (RF)
#Libraries
library(ranger)             # Random forest
library(tidyverse)
library(hydroGOF)           # rmse()
library(caret)              # Train
library(StratifiedMedicine) #Plot importance
plot()
#functions
source("E:\\zuomian\\毕业论文\\论文修改\\盆地污染物估计\\各种代码\\plot_importance.R")
source("E:\\zuomian\\毕业论文\\论文修改\\盆地污染物估计\\各种代码\\nrmse_func.R")
df_clean <- data.frame(
  status = df$out,
  time = df$`Survival months`,
  Age = df$agegroup,
  Histology = df$histology,
  NumberT = df$`Total number of in situ/malignant tumors for patient`,
  Continuum = df$Continuum,
  Income = df$income,
  MaritalStatus = df$`Marital status at diagnosis`,
  NodesExamined = df$`Regional nodes examined (1988+)`,
  NodesPositive = df$`Regional nodes positive (1988+)`,
  Chemotherapy = df$`Chemotherapy recode (yes, no/unk)`,
  Sex = df$Sex,
  Race = df$`Race recode (White, Black, Other)`,
  Radiation=df$radiation
 # time=df$time
)
df_clean$status <- as.factor(df_clean$status)
df_clean <- as.data.frame(lapply(df_clean, function(x) {
  if(is.character(x)) factor(x) else x
}))
# 划分训练集和测试集
set.seed(123)
table(df_clean$status)
df_clean$status <- as.numeric(df_clean$status)-1
# 构建随机森林训练模型
rf_model <- rfsrc(Surv(time, status) ~ ., data = df_clean,nsplit = 10, 
                  na.action = "na.impute", tree.err = TRUE,  
                  importance = TRUE, ntree = 300)
get.cindex(time = df_clean$time, censoring = df_clean$status, predicted = rf_model$predicted.oob)
srv.smp.o <- subsample(rf_model, B = 100, performance = TRUE)
plot(srv.smp.o)
# 生成生存概率预测
preds <- predict.rfsrc(rf_model, newdata = df_clean)
predicted_survival_5_years <- as.vector(preds$survival[,36])
library(pROC)
# 创建ROC对象
roc_obj <- roc(df_clean$status, predicted_survival_5_years,xlim =c(0,1))
# 绘制ROC曲线
windowsFonts(roman = windowsFont("TT TIMES"))
par(family = "roman")
plot(roc_obj,
     print.thres=TRUE, # 输出cutoff值
     print.thres.cex=0.8, # 设置cutoff值字体大小
     col= "#377eb8", # 曲线颜色
     print.thres.col="#e41a1c", # cutoff值字体的颜色
     identity.col="#4daf4a", # 对角线颜色
     identity.lty=2, identity.lwd=2, # 对角线样式
     main = "ROC Curve for Survival at 5 Years",
     legacy.axes=TRUE)

# 计算AUC值和置信区间
auc_value <- auc(roc_obj)
auc_ci <- ci(roc_obj)

# 在图上标记最佳灵敏度和特异性点
best_coords <- coords(roc_obj, "best", ret=c("sensitivity", "specificity", "threshold"))
#points(best_coords[1, "threshold"], best_coords[1, "sensitivity"], pch=19, col="#ff7f00") # 使用橙色标记点
# 添加图例
legend("bottomright",
       legend=paste0("AUC = ", round(auc_value, 3), 
                     " [", round(auc_ci[1], 3), "-", round(auc_ci[2], 3), "]",
                     "\nBest (Sensitivity, Specificity) = (", round(best_coords[1, "sensitivity"], 3), ", ", round(best_coords[1, "specificity"], 3), ")"),
       box.lty=0,
       col="#000000", # 文本颜色
       cex=0.8, # 字体大小
       bg="lightgrey") # 图例背景颜色


