matrix_data <- matrix(data=results_hglasso$`Case4_p=100`$Hub_accuracy_node,nrow = length(seq(0.02, 0.3, 0.02)))

matrix_df <- melt(matrix_data)
 # 绘制热图
 ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
       geom_tile() +
       scale_fill_gradient(low = "blue", high = "red") +
       theme_minimal() +
   labs(title = "",
        x = paste0("Combinations of ",expression(lambda_2)," and ",expression(lambda_3)),
        y = expression(lambda_1),
        color = "Methods") +
   theme(
     plot.title = element_text(hjust = 0.5, size = 36), # 放大标题
     axis.title.x = element_text(size = 26),             # 放大x轴标签
     axis.title.y = element_text(size = 30),             # 放大y轴标签
     axis.text.x = element_text(size = 24),              # 放大x轴数值
     axis.text.y = element_text(size = 24),              # 放大y轴数值
     legend.title = element_text(size = 30),             # 放大图例标题
     legend.text = element_text(size = 30),
     legend.spacing = unit(3, "cm"),                 # 图例间距
     legend.key.height = unit(2.5, "cm") # 放大图例内容
   )
 