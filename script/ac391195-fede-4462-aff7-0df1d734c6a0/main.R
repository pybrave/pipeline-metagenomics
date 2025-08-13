library(tidyverse)
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)

print(args)
params_path <- args[1]
output_path <- args[2]


data <- fromJSON(params_path)
control <- data$control |>
    mutate(select_group= data$groups_name$control)
treatment <- data$treatment |>
    mutate(select_group= data$groups_name$treatment)

list_path <- rbind(control, treatment)

metadata <- list_path[c("sample_name","select_group")]


df_list <- apply(list_path,1, function(x){
    term_path <- x[[data$term]]
    sample_name <- x[["sample_name"]]
    df <-  read_tsv(term_path)
    colnames(df) <- c("term",sample_name)
    df <- df |> filter(!grepl("\\|", term))
    df
    # print()
})

merger_df <- reduce(df_list,function(x,y){merge(x,y,by="term",all=T)} ) |>
    column_to_rownames("term") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))


merger_df <- merger_df[metadata$sample_name]

result <- data.frame(
  term = rownames(merger_df),
  p_value = apply(merger_df, 1, function(vals) {
    wilcox.test(as.numeric(vals) ~ metadata$select_group)$p.value
  })
)
pseudo_count <- 1e-6
log2fc <- apply(merger_df, 1, function(vals) {
  mean1 <- mean(vals[metadata$select_group == data$groups_name$treatment], na.rm = TRUE)
  mean2 <- mean(vals[metadata$select_group == data$groups_name$control], na.rm = TRUE)
    # print(c(mean1,mean2,log2(mean1+pseudo_count / mean2+pseudo_count)))
  log2(mean1+pseudo_count / mean2+pseudo_count)
})
result$log2fc <- log2fc

result$FDR <- p.adjust(result$p_value, method = "BH")
result <- result[order(result$FDR), ]
filename <- paste0(unique( metadata$select_group),collapse = "_vs_")
result  |>
    left_join(rownames_to_column(merger_df,"term"),by="term") |>
    write_tsv(paste0(output_path,"/",filename,".tsv"))