
# while (TRUE) {
#   Sys.sleep(1000)  # 每次休眠 1000 秒，然后继续循环
# }

library(Maaslin2)
library(jsonlite)
library(tidyverse)
library(pheatmap)
# library(ggrepel)
library(ggrepel)   # 用于防止标签重叠
library(lefser)
library(ggpubr)


args <- commandArgs(trailingOnly = TRUE)
print(args)

params_path <- args[1]
output_path <- args[2]
if(F){
  params_path <- "params.json"
  output_path <- "output"
}


data <- fromJSON(params_path)
control <- data$control |>
  mutate(select_group= data$groups_name$control)
treatment <- data$treatment |>
  mutate(select_group= data$groups_name$treatment)

list_path <- rbind(control, treatment)

metadata <- list_path[c("sample_name","select_group")]

rank <- data$rank

duplicated(list_path$sample_name)


read_abundabce <- function(path){
  df <-  read_tsv(path,comment = "#",col_names =F)
  colnames(df) <- c("clade_name","NCBI_tax_id","abundance","additional_species")
  df <- select(df,c("clade_name",all_of("abundance")))
  df
}
parse_metaphlan <- function(df,sample_name,rank) {
  df %>%
    mutate(clade_name = as.character(clade_name)) %>%
    separate_wider_delim(
      clade_name,
      delim = "|",
      names = c("KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES","SGB"),
      too_few = "align_start"
    )|>
    mutate(row_rank=case_when(!is.na(SGB)~'SGB',
                          !is.na(SPECIES)~'SPECIES',
                          !is.na(GENUS)~'GENUS',
                          !is.na(FAMILY)~'FAMILY',
                          !is.na(ORDER)~'ORDER',
                          !is.na(CLASS)~'CLASS',
                          !is.na(PHYLUM)~'PHYLUM',
                          !is.na(KINGDOM)~'KINGDOM'))|>
    mutate(sample_name= sample_name) |> 
    filter(row_rank==rank) |>
    select(sample_name,taxonomy=all_of(rank),abundance) 
    
}

# sample_name <- "OCC8"
# df <- read_abundabce("/ssd1/wy/workspace2/nextflow_workspace/289364b1-295c-4710-833e-d68ec7c8918e/131f8806-35e3-4d7c-b234-f14a2119aaa7/2c88b345-822f-4285-9222-a18b9c3daa8b/output/metaphlan/OCC8/OCC8_profile.txt")
# a <- parse_metaphlan(df,"aa","SPECIES") |>
#   filter(!grepl("GGB|GBS",taxonomy)) |>
#   mutate(abundance= abundance/sum(abundance)*100)
# sum(a$abundance)

df_list <- apply(list_path,1, function(x){
  profile_path <- x[["profile"]]
  sample_name <- x[["sample_name"]]
  # df <-  read_tsv(term_path,comment = "#")
  # colnames(df) <- c("clade_name",sample_name)
  df <- read_abundabce(profile_path)
  df <- parse_metaphlan(df,sample_name,rank) 
  if(data$filter_unknown_taxonomy){
    df <- filter(df,!grepl("GGB|GBS",taxonomy)) |>
      mutate(abundance= abundance/sum(abundance)*100)
  }
  # df <- df |> filter(!grepl("\\|", term))
  df
})

df_long <- bind_rows(df_list)

merged_df <- df_long %>%
  pivot_wider(names_from = sample_name, values_from = abundance) |>
  mutate(across(where(is.numeric), ~replace_na(., 0))) |>
  column_to_rownames("taxonomy") 


# merged_df
metadata <- metadata |>column_to_rownames("sample_name")
# if your inputs are counts, then you can use NEGBIN and ZINB, 
# whereas, for non-count (e.g. percentage, CPM, or relative abundance) input, you can use LM and CPLM.

fit_data = Maaslin2(input_data     = t(merged_df) |> as.data.frame(), 
                    input_metadata = metadata, 
                    plot_scatter =data$plot_scatter,
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "output", 
                    fixed_effects  = c("select_group"),
                    reference      = c(paste0("select_group,",data$groups_name$control)))  

all_results <-fit_data$results  #read_tsv("output/all_results.tsv")
# read_tsv("output/significant_results.tsv") 

title <- paste0(c(data$groups_name$treatment,data$groups_name$control),collapse = " vs ")



se <- SummarizedExperiment(
  assays = list(count=as.matrix(merged_df)),
  colData =metadata
)

set.seed(1234)
setn_ra <- relativeAb(se)
b <- assay(setn_ra)
res1 <- lefser(setn_ra, # relative abundance only with terminal nodes
               kruskal.threshold=1,
               lda.threshold=0,
               classCol = "select_group")
# dim(res1)
# dim(all_results)
df_res <-  all_results |>
  left_join(dplyr::rename(as.data.frame(res1),"feature"=features),by="feature")

df_res |> write_tsv(file = paste0(output_path,"/",str_replace_all(title," ","_"),".score.tsv"))



# 转换 qval 为 -log10
sig_thresh <-  data$sig_thresh
effect_cutoff <- data$effect_cutoff
label_size <- data$label_size

all_results <- all_results %>%
  mutate(sig_value = .data[[data$sig_type]]) |>
  mutate(direction = factor(ifelse(sig_value  < sig_thresh & abs(coef) > effect_cutoff ,
                                   ifelse(coef>0,"Up","Down"),"NS"),
                            levels = c("Up","Down","NS") 
  )) 
counts <- table(all_results$direction)
labels <- c(
  Down = paste0("Down (", counts["Down"], ")"),
  NS   = paste0("NS (", counts["NS"], ")"),
  Up   = paste0("Up (", counts["Up"], ")")
)

ggplot(all_results, aes(x= coef, 
                            y = -log10(sig_value), 
                            colour=direction)) +
  geom_point(alpha=0.9, size=3.5)+
  scale_color_manual(values=c("Down" = "#3B4992FF", "NS"="#d2dae2","Up"="#EE0000FF"),  labels = labels)+
  geom_vline(xintercept=c(-effect_cutoff,effect_cutoff),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(sig_thresh),lty=4,col="black",lwd=0.8)+
  ggtitle(paste0("volcano plot of ",title)) +
  labs(x="Effect Size (Coefficient)", y="-log10(q-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(color = 'black',size = 10, family = 'sans', face = 'plain'),
        axis.title = element_text(color = 'black',size = 15, family = 'sans', face = 'plain'),
        legend.position="right", 
        legend.title = element_blank()
  ) + 
  geom_text_repel(data=all_results %>% filter(direction!="NS") |> head(label_size),aes(label = feature), 
                  colour = "black", size = 4)


ggsave(filename = paste0(output_path,"/",str_replace_all(title," ","_"),".pdf"))

feature_list <- all_results |> arrange(pval) |>pull(feature)
df_long |>
  filter(taxonomy == feature_list[1]) |>
  inner_join(rownames_to_column(metadata,"sample_name"),by="sample_name") |>
  mutate(select_group =factor(select_group, 
                              levels = c(data$groups_name$treatment,data$groups_name$control)) ) |>
  ggplot(aes(x=select_group, y=abundance )) +
  geom_boxplot(fill = "skyblue", color = "black") +
  geom_point()+
  theme_bw()+
  ggtitle(paste0(all_results[all_results$feature==feature_list[1],c("feature")],
                 "\n",title," coef: ",
                 all_results[all_results$feature==feature_list[1],c("coef")]))
ggsave(filename = paste0(output_path,"/",str_replace_all(title," ","_"),".boxplot.pdf"))

ann_colors <- list(
  group = setNames(
    c("steelblue", "tomato"),  # 颜色向量
    c( data$groups_name$control, data$groups_name$treatment)  # 用变量的值作为名字
  )
)
dev.off()
pdf(file =paste0(output_path,"/",str_replace_all(title," ","_"),".heatmap.pdf") , width = data$heatmap_width,height =8)

log2(merged_df[head(feature_list, n=30),]+1) |>
  pheatmap(scale = "row", 
           annotation_colors = ann_colors,
           cluster_cols = data$heatmap_cluster_cols,
           show_colnames = data$heatmap_show_colnames,
           color =  colorRampPalette(c("darkred", "#FFFFFF","darkblue"))(255),
           annotation_col=metadata |> dplyr::rename(group=select_group)
           )
dev.off()




sig_feature <- df_res |>
  arrange(qval) |>
  head(n=data$boxplot_num) |>
  mutate(sig_value = .data[[data$sig_type]]) |>
  mutate(p.signif = case_when(
    sig_value < 0.001 ~ "***",
    sig_value < 0.01  ~ "**",
    sig_value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))
box_data <-  df_long |>
  mutate(abundance=log2(abundance+data$boxplot_pseudo_count))
  
# log2(t(merged_df)[head(feature_list, n=30),]+1) 
p <- box_data |>
  filter(taxonomy %in% pull(sig_feature,"feature")) |>
  inner_join(rownames_to_column(metadata,"sample_name"),by="sample_name") |>
  ggplot( aes(x=taxonomy, y=abundance, fill=select_group)) +
  geom_boxplot(position=position_dodge(0.8), width=0.6, outlier.shape=NA) + # 去掉离群点，箱子更美观
  geom_jitter(aes(color=select_group), 
              position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8), 
              alpha=0.5, size=1.5) +  # 添加散点，更直观
  labs(title=paste0("Boxplot of ",title), 
       x="", 
       y="log2(abundance)") +
  scale_fill_brewer(palette="Set2") +   # 柔和调色板
  scale_color_brewer(palette="Set2") +  # 散点与箱子颜色一致
  theme_bw(base_size=12) +
  theme(
    axis.text.x = element_text(angle=80, vjust=1, hjust=1, size=10, face="italic"), 
    axis.title.x = element_text(size=12, face="bold"),  
    axis.title.y = element_text(size=12, face="bold"),  
    plot.title = element_text(hjust=0.5, size=14, face="bold"),  
    panel.grid.major.x = element_blank(),   # 去掉竖网格线
    panel.grid.minor = element_blank(),  
    legend.title = element_blank(),        # 去掉图例标题
    legend.position = "top"                # 图例放上方
  )

if (data$boxplot_sig_label =="Maaslin2"){
  p=p+geom_text(
    data = sig_feature,
    aes(x = feature, y = max(box_data$abundance), 
        label = p.signif),
    inherit.aes = FALSE
  ) 
}else{
  p=p+stat_compare_means(
    aes(group=select_group),
    method="wilcox.test",   # 或 "t.test"
    label="p.signif",       # 显示星号，也可用 "p.format" 显示具体 p 值
    hide.ns=TRUE           # 不显示 ns
  )
}
p


ggsave(filename = paste0(output_path,"/",str_replace_all(title," ","_"),"multigroup.boxplot.pdf"),plot = p,width = 15,height = 8)

#   mutate(log10q = -log10(qval),
#          signif = qval < 0.0001 & abs(coef)>5)  # 标记显著
# 
# # 火山图
# ggplot(all_results, aes(x = coef, y = log10q, color = signif)) +
#   geom_point(size = 3, alpha = 0.8) +
#   scale_color_manual(values = c("grey", "red")) +
#   theme_minimal() +
#   labs(x = "Effect Size (Coefficient)", y = "-log10(q-value)",
#        title = "MaAsLin 2 Volcano Plot") +
#   geom_text(data = all_results %>% filter(signif),
#             aes(label = feature),
#             vjust = 1.5, hjust = 0.5, size = 3)
# 
# volcano <- function(res,title){
#   pval_cutoff = 0.05
#   logfc_cutoff = 1
#   deg_res_volcano <- as.data.frame(res) |>
#     # rownames_to_column("symbol") |> 
#     mutate(direction = factor(ifelse(pvalue < pval_cutoff & abs(log2FoldChange)>logfc_cutoff,
#                                      ifelse(log2FoldChange>0,"Up","Down"),"NS"),
#                               levels = c("Up","Down","NS") 
#     )) |>
#     mutate(label =ifelse(pvalue < pval_cutoff & abs(log2FoldChange)> logfc_cutoff, Symbol,"") ) 
#   
#   ggplot(deg_res_volcano, aes(x= log2FoldChange, 
#                               y = -log10(pvalue), 
#                               colour=direction)) +
#     geom_point(alpha=0.9, size=3.5)+
#     scale_color_manual(values=c("Down" = "#3B4992FF", "NS"="#d2dae2","Up"="#EE0000FF"))+
#     geom_vline(xintercept=c(-logfc_cutoff,logfc_cutoff),lty=4,col="black",lwd=0.8) +
#     geom_hline(yintercept = -log10(pval_cutoff),lty=4,col="black",lwd=0.8)+
#     ggtitle(paste0("volcano plot of ", title)) +
#     labs(x="log2(fold change)", y="-log10 (p-value)")+
#     theme_bw()+
#     theme(plot.title = element_text(hjust = 0.5), 
#           axis.text = element_text(color = 'black',size = 10, family = 'sans', face = 'plain'),
#           axis.title = element_text(color = 'black',size = 15, family = 'sans', face = 'plain'),
#           legend.position="right", 
#           legend.title = element_blank()
#     )+ geom_text_repel(data=deg_res_volcano,aes(label = label), 
#                        colour = "black", size = 4)
# }





# input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
# input_data
# 
# input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
# input_metadata
# 
# df_input_data = read.table(file             = input_data,
#                            header           = TRUE,
#                            sep              = "\t", 
#                            row.names        = 1,
#                            stringsAsFactors = FALSE)
# df_input_data[1:5, 1:5]
# 
# df_input_metadata = read.table(file             = input_metadata, 
#                                header           = TRUE, 
#                                sep              = "\t", 
#                                row.names        = 1,
#                                stringsAsFactors = FALSE)
# df_input_metadata[1:5, ]
# 
# 
# fit_data = Maaslin2(input_data     = input_data, 
#                     input_metadata = input_metadata, 
#                     min_prevalence = 0,
#                     normalization  = "NONE",
#                     output         = "demo_output", 
#                     fixed_effects  = c("diagnosis"),
#                     reference      = c("diagnosis,nonIBD"))

# 
# df_long <- map_dfr(seq_len(nrow(list_path)), function(i) {
#   profile_path <- list_path$profile[i]
#   sample_name  <- list_path$sample_name[i]
#   
#   df <- read_abundabce(profile_path, sample_name) %>%
#     parse_metaphlan(sample_name, rank) %>%
#     mutate(sample = sample_name) %>%
#     rename(abundance = all_of(sample_name))
#   
#   df %>% select(sample, taxonomy, abundance)
# })
