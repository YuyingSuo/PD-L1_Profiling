library(ggplot2)
library(ggpubr)

ggplot_theme=function(){
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )
}



#@@@@@@@@@@@@@@@@ expression

## NAT & T
draw_expr_plot <- function(data) {
  n=data.frame(t(data[,grep("_N", colnames(data))]))
  t=data.frame(t(data[,-grep("_N", colnames(data))]))
  pval <- wilcox.test(as.numeric(n$CD274),as.numeric(t$CD274))$p.value
  pval_label <- paste0("wilcox, p = ", formatC(pval, format = "e", digits = 2))
  
  n$sample='NAT'
  t$sample='T'
  tmp2=rbind(n,t)
  y_max <- max(tmp2$CD274, na.rm = TRUE)
  p <-ggplot(tmp2, aes(x=sample, y=CD274, color=sample)) +
    scale_color_manual(values = c("#55A0FB","#FF8080"))+
    geom_boxplot(width = 0.7,outlier.size = 1)+
    annotate("text", x = 1.5, y = y_max * 1.05, label = pval_label, size = 5, color = "black") +
    labs(title = "PD-L1 expression")+
    ggplot_theme()

  return(p)
}


## only T
draw_expr_plot_onlyT <- function(data){
  tmp2=t(data)
  p=ggplot(tmp2, aes(x="", y=CD274, color="#FF8080")) +
    geom_boxplot(width = 0.5,outlier.size = 1)+
    theme(axis.text.x =element_text(angle=60,color = "black",hjust = 1), 
          axis.text.y=element_text(color = "black"))+
    labs(x='T',title = "PD-L1 expression")+
    ggplot_theme()
  return(p)
}





#@@@@@@@@@@@@@@@@ DEG
draw_deg_plot <- function(data,gene){
  deg_sig_gene=data[which(rownames(data) == gene),]
  
  dt=deg_sig_gene[, -((ncol(deg_sig_gene)-3):ncol(deg_sig_gene))]
  dt2=data.frame(t(dt))
  try(colnames(dt2) <- 'value', silent = TRUE)
  dt2$stage=sub(".*_", "", rownames(dt2))
  #print(head(dt2))
   p = NULL
  try({
    pval <- wilcox.test(as.numeric(dt2[which(dt2$stage=='low'),1]),as.numeric(dt2[which(dt2$stage=='high'),1]))$p.value
    pval_label <- paste0("wilcox, p = ", formatC(pval, format = "e", digits = 2))
    y_max <- max(dt2[,1], na.rm = TRUE)
    
    p = ggplot(dt2, aes(x=stage, y=value, color=stage)) +
      scale_color_manual(
        values = c("low" = "#55A0FB", "high" = "#FF8080"),
        labels = c("low" = "PD-L1_Low", "high" = "PD-L1_High")
      )+
      scale_x_discrete( labels = c("low" = "PD-L1_Low", "high" = "PD-L1_High")) +
      geom_boxplot(width = 0.7,outlier.size = 1)+
      labs(y=paste0(gene,"_exprssion"))+
      annotate("text", x = 1.5, y = y_max * 1.05, label = pval_label, size = 5, color = "black")+
      ggplot_theme()
  }, silent = TRUE)
  
  return(p)
}





#@@@@@@@@@@@@@@@@ Cor
draw_cor_plot <- function(data,gene){
  cor_select=data[,-((ncol(data)-3):ncol(data))]
  #print(head(cor_select))
  cor_tmp=data.frame(t(cor_select))

  p = ggscatter(cor_tmp,x=names(cor_tmp)[2],y=names(cor_tmp)[1],add = "reg.line",   size=1 ,              
            conf.int = TRUE,                    
            add.params = list(color = "black",  
                              fill = "#F9D0A0"),
            ggtheme = ggplot_theme()) +
    stat_cor(method = "spearman")
  return(p)
}





#@@@@@@@@@@@@@@@@ path
draw_path_plot <- function(data){
  data=data.frame(t(data))
  n=data.frame(t(data[,grep("_N", colnames(data))]))
  t=data.frame(t(data[,-grep("_N", colnames(data))]))
  pval <- wilcox.test(as.numeric(n$ssgsea_PD_L1_path),as.numeric(t$ssgsea_PD_L1_path))$p.value
  pval_label <- paste0("wilcox, p = ", formatC(pval, format = "e", digits = 2))
  
  n$sample='NAT'
  t$sample='T'
  tmp2=rbind(n,t)
  
  y_max <- max(tmp2$ssgsea_PD_L1_path, na.rm = TRUE)
  p=ggplot(tmp2, aes(x=sample, y=ssgsea_PD_L1_path, color=sample)) +
    scale_color_manual(values = c("#55A0FB","#FF8080"))+
    geom_boxplot(width = 0.7,outlier.size = 1)+
    annotate("text", x = 1.5, y = y_max * 1.05, label = pval_label, size = 5, color = "black") +
    labs(title = "PD-L1 pathway ssgsea scores")+
    ggplot_theme()
  return(p)
}

draw_path_plot_onlyT <- function(data){
  p = ggplot(data, aes(x="", y=ssgsea_PD_L1_path, color="#FF8080")) +
    geom_boxplot(width = 0.5,outlier.size = 1)+
    theme(axis.text.x =element_text(angle=60,color = "black",hjust = 1), 
          axis.text.y=element_text(color = "black"))+
    labs(x='T',title = "PD-L1 pathway ssgsea scores")+
    ggplot_theme()
  return(p)
}






#@@@@@@@@@@@@@@@@ pho
draw_pho_plot <- function(data,site){
  pho=data.frame(t(data[,-((ncol(data)-3):ncol(data))]))
  p = NULL
  try({
    colnames(pho)="value"
    pho$stage=sub(".*_", "", rownames(pho))
    
    library(ggplot2)
    y_max <- max(pho$value, na.rm = TRUE)
    pval <- data[site,"p"]
    pval_label <- paste0("wilcox, p = ", formatC(pval, format = "e", digits = 2))
    
    p <- ggplot(pho, aes(x=stage, y=value, color=stage)) +
      scale_color_manual(
        values = c("Inactive" = "#55A0FB", "Active" = "#FF8080"),
        labels = c("Inactive" = "PD-L1_Path Inactive", "Active" = "PD-L1_Path Active")
      )+
      geom_boxplot(width = 0.7,outlier.size = 1)+
      annotate("text", x = 1.5, y = y_max * 1.05, label = pval_label, size = 5, color = "black") +
      labs(title = "PhosphoSites Abundance")+
      ggplot_theme()
  }, silent = TRUE)
  
  return(p)
}