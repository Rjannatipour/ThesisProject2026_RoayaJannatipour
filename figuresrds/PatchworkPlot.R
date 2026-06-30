#Patchwork plot
install.packages("ggVennDiagram")
install.packages("patchwork")

library(ggVennDiagram)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(flextable)
library(readr)
library(officer)




#Males

#Univariate Plots#
#Random Forest
rfmaleunivariateplot <- readRDS("figuresrds/rfmaleunivariateplot.rds") + ggtitle("Univariate + RF ") + geom_point(color = "blue", alpha=0.1)

#Elastic Net
elasticnetmaleunivariateplot <- readRDS("figuresrds/elasticnetmaleunivariateplot.rds") + ggtitle("Univariate + Elastic Net ") + geom_point(color = "blue", alpha = 0.1)

#RISE#
#Random Forest

rfmaleriseplot <- readRDS("figuresrds/rfmaleriseplot.rds") + ggtitle("RISE + RF ") + geom_point(color = "blue", alpha = 0.1)

#Elastic Net
elasticnetmaleriseplot <- readRDS("figuresrds/elasticnetmaleriseplot.rds")+ ggtitle("RISE + Elastic Net ") + geom_point(color = "blue", alpha = 0.1)




#Females

#Univariate Plots#
#Random Forest
rffemaleunivariateplot <- readRDS("figuresrds/rfemaleunivariateplot.rds") + ggtitle("Univariate + RF") 
rffemaleunivariateplot$layers[[1]] <- 
  geom_point(color = "pink", alpha = 0.6)

#Elastic Net
elasticnetfemaleunivariateplot <- readRDS("figuresrds/elasticnetfemaleunivariateplot.rds")+ ggtitle("Univariate + Elastic Net") 
elasticnetfemaleunivariateplot$layers[[1]] <- 
  geom_point(color = "pink", alpha = 0.6)

#RISE#
#Random Forest
rffemaleriseplot <- readRDS("figuresrds/rffemaleriseplot.rds")+ ggtitle("RISE + RF") 
rffemaleriseplot$layers[[1]] <- 
  geom_point(color = "pink", alpha = 0.6)

#Elastic Net#
elasticnetfemaleriseplot <- readRDS("figuresrds/elasticnetfemaleriseplot.rds")+ ggtitle("RISE + Elastic Net") 
elasticnetfemaleriseplot$layers[[1]] <- 
  geom_point(color = "pink", alpha = 0.6)






# Combine
combined_plot_male <- (rfmaleunivariateplot | rfmaleriseplot) /
  (elasticnetmaleunivariateplot        | elasticnetmaleriseplot) + plot_annotation(theme = theme(
    plot.title = element_text(size = 18)
  )) + plot_annotation(title = 'Male Participants - Prediction Methods')

ggsave(
  "figurespng/combined_plot_male.png",
  combined_plot_male,
  width = 6,
  height = 6,
  dpi = 300
)

combined_plot_female <- (rffemaleunivariateplot | rffemaleriseplot) /
  (elasticnetfemaleunivariateplot  | elasticnetfemaleriseplot)+ plot_annotation(theme = theme(
    plot.title = element_text(size = 18)
  )) + plot_annotation(title = 'Female Participants - Prediction Methods')


ggsave(
  "figurespng/combined_plot_female.png",
  combined_plot_female,
  width = 6,
  height = 6,
  dpi = 300
)



# Venn Diagram - Univariate Screening

justimprovegenesfemale <- read.csv("FemaleUnivariateScreening/justimprovegenesfemale.csv")

justimprovegenesfemale <- justimprovegenesfemale %>% pull(Genes)

justimprovegenesmale <- read.csv("MaleUnivariateScreening/justimprovegenesmale.csv")

justimprovegenesmale <- justimprovegenesmale %>% pull(Genes)


univenn <- list(
  Female = justimprovegenesfemale,
  Male   = justimprovegenesmale
)

univenndiagram <- ggVennDiagram(univenn) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void() + coord_flip()

ggsave(
  "figurespng/univenndiagram.png",
  univenndiagram,
  width = 8,
  height = 6,
  dpi = 300
)

#RISE Venn Diagram
risesigfemale <- readRDS("FemaleRISEPrediction/risefemalesig.rds")

risesigmale <- readRDS("MaleRISEPrediction/risemalesig.rds")

risevenn <- list(
  Female = risesigfemale,
  Male   = risesigmale
)

risevenndiagram <- ggVennDiagram(risevenn) +
  scale_fill_gradient(
    
    low = "white",   # white
    high = "#08306b"   # deep navy
  ) +
  theme_void() + coord_flip()

ggsave(
  "figurespng/risevenndiagram.png",
  risevenndiagram,
  width = 8,
  height = 6,
  dpi = 300
)

# Common RISE genes table

common_genes <- intersect(risesigfemale, risesigmale)
half <- 10
common_df <- data.frame(
  `Gene (1-10)`  = common_genes[1:half],
  `Gene (11-20)` = common_genes[(half+1):length(common_genes)],
  check.names = FALSE
)

commonrisegenes <- flextable(common_df) %>%
  set_caption(caption = "Common genes after RISE screening in  Males and Females") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  vline(j = 1, border = fp_border(width = 0.5)) %>%
  autofit()


save_as_html(commonrisegenes,
             path = "figurespng/commonrisegenes.html")

#Presenting Common Univariate Screening Genes:

commonuni <- intersect(justimprovegenesmale,justimprovegenesfemale)

commonunitop20 <- commonuni[1:20]

common_df_uni <- data.frame(
  `Gene (1-10)`  = commonunitop20[1:10],
  `Gene (11-20)` = commonunitop20[11:20],
  check.names = FALSE
)

commonunigenes <- flextable(common_df_uni) %>%
  set_caption(caption = "Common Genes after Univariate Screening in Males and Females") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  vline(j = 1, border = fp_border(width = 0.5)) %>%
  autofit()


save_as_html(commonunigenes,
             path = "figurespng/commonunigenes.html")

# Top 10 Performing Genes Univariate Model - Males

top10malegenes <- read_csv("MaleUnivariateScreening/top10male.csv")

top10maletable <- flextable(top10malegenes) %>%
  set_caption(caption = "Top 10 Performing Genes in Univariate Gene Model - Males") %>%
  set_header_labels(Genes = "Gene", `Diff_RMSE` = "Difference sRMSE") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  autofit()

save_as_html(top10maletable,
             path = "MaleUnivariateScreening/top10maletable.html")

# Top 10 Performing Genes Univariate Model - Females

top10femalegenes <- read_csv("FemaleUnivariateScreening/top10female.csv")

top10femaletable <- flextable(top10femalegenes) %>%
  set_caption(caption = "Top 10 Performing Genes in Univariate Gene Model - Females") %>%
  set_header_labels(Genes = "Gene", `Diff_RMSE` = "Difference sRMSE") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  autofit()

save_as_html(top10femaletable,
             path = "FemaleUnivariateScreening/top10femaletable.html")





