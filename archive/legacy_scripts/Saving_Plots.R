load("F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/Results1and2.RData")
load("F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/Results3_0.5.RData")

load("F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/Results_4_hub=0.2.RData")

`Case1_p=100`



# 安装并加载ggplot2包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

datas <- list(`Case1_p=100`,`Case1_p=200`,`Case1_p=500`,`Case2_p=100`,`Case2_p=200`,`Case2_p=500`)
vari <- c("BIC","TPR","FPR","Sparsity_ratio","Fnorm","Time")
name <- c("_1_100","_1_200","_1_500","_2_100","_2_200","_2_500")


for (i in 1:6) {
  
  for (j in vari) {
    
    
    ydata <- datas[[i]]
    
    
    # 示例数据
    x <- seq(0.1,2,0.05)  # 横坐标向量
    y1 <- ydata[[j]]$glasso[-1]  # 第一个向量
    y2 <- ydata[[j]]$dpglasso[-1]  # 第二个向量
    y3 <- ydata[[j]]$QUIC[-1]  # 第三个向量
    
    # 整理数据为长格式
    data <- data.frame(
      x = rep(x, 3),
      y = c(y1, y2, y3),
      group = rep(c("glasso", "dpglasso", "QUIC"), each = length(x))
    )
    
    
    
    
    # 绘图
    ggplot(data, aes(x = x, y = y, color = group)) +
      geom_line(linewidth = 1) +  # 折线
      geom_point(size = 2) + # 数据点
      labs(title = j,
           x = expression(lambda),
           y = j,
           color = "Methods") + 
      theme_minimal()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 36), # 放大标题
        axis.title.x = element_text(size = 36),             # 放大x轴标签
        axis.title.y = element_text(size = 30),             # 放大y轴标签
        axis.text.x = element_text(size = 24),              # 放大x轴数值
        axis.text.y = element_text(size = 24),              # 放大y轴数值
        legend.title = element_text(size = 36),             # 放大图例标题
        legend.text = element_text(size = 36),
        legend.spacing = unit(2, "cm"),                 # 图例间距
        legend.key.height = unit(2, "cm") # 放大图例内容
      )
    
    save_path <- "F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/plots"
    
    
    
    file_name <- file.path(save_path, paste0(j, name[i], ".png"))
    
    ggsave(file_name, width = 8, height = 6, dpi = 300)
    
    
  }
  
}

load("F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/Results1and2.RData")
load("F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/Results3_0.5.RData")




# 安装并加载ggplot2包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

datas <- list(`Case3_p=100_Sp=0.5`,`Case3_p=200_Sp=0.5`,`Case3_p=500_Sp=0.5`,`Case3_p=100_Sp=0.99`,`Case3_p=200_Sp=0.99`,`Case3_p=500_Sp=0.99`)
vari <- c("BIC","TPR","FPR","Sparsity_ratio","Fnorm","Time")
name <- c("_3_100_0.5","_3_200_0.5","_3_500_0.5","_3_100_0.99","_3_200_0.99","_3_500_0.99")


for (i in 1:6) {
  
  for (j in vari) {
    
    
    ydata <- datas[[i]]
    
    
    # 示例数据
    x <- seq(0.1,2,0.05)  # 横坐标向量
    y1 <- ydata[[j]]$glasso  # 第一个向量
    y2 <- ydata[[j]]$dpglasso  # 第二个向量
    y3 <- ydata[[j]]$QUIC  # 第三个向量
    
    # 整理数据为长格式
    data <- data.frame(
      x = rep(x, 3),
      y = c(y1, y2, y3),
      group = rep(c("glasso", "dpglasso", "QUIC"), each = length(x))
    )
    
    
    
    
    # 绘图
    ggplot(data, aes(x = x, y = y, color = group)) +
      geom_line(linewidth = 1) +  # 折线
      geom_point(size = 2) + # 数据点
      labs(title = j,
           x = expression(lambda),
           y = j,
           color = "Methods") + 
      theme_minimal()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 36), # 放大标题
        axis.title.x = element_text(size = 36),             # 放大x轴标签
        axis.title.y = element_text(size = 30),             # 放大y轴标签
        axis.text.x = element_text(size = 24),              # 放大x轴数值
        axis.text.y = element_text(size = 24),              # 放大y轴数值
        legend.title = element_text(size = 36),             # 放大图例标题
        legend.text = element_text(size = 36),
        legend.spacing = unit(2, "cm"),                 # 图例间距
        legend.key.height = unit(2, "cm") # 放大图例内容
      )
    
    save_path <- "F:/UCD Courses/2024 Fall/STA 250/Final Project/GLASSO/results/plots"
    
    
    
    file_name <- file.path(save_path, paste0(j, name[i], ".png"))
    
    ggsave(file_name, width = 8, height = 6, dpi = 300)
    
    
  }
  
}

`Case2_p=500`
`Case3_p=500_Sp=0.99`
