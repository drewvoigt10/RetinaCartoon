#' Generate data for the retina cartoon
#'
#' \code{generate_cartoon_data} returns a list with two elements. The first element is a
#' heatmap-cartoon genenerated by the grImport package. The second is a legend that related
#' expression to color. These elements can then be combined by the plot_cartoon function
#' to display a pretty cartoon plot
#'
#' @param gene Character string of the gene name for expression visualization. If the requested
#' gene is not found, it was most likely not annotated in the given dataset. Try adding the suffix
#' '-GRCh38' if you believe the gene may have been annotated in the more recent genome build.
#'
#' @param dataset Character string of the dataset name. Options are: "godzilla", "retina_fov_vs_periph",
#' "retina_AIR_vs_control", "RPE_choroid_unselected", "RPE_choroid_CD31_selected", and
#' "CD31_choroid_infant_adult". Godzilla (a merged object of all retina+rpe+choroid
#' single-cell datasets published by 09/01/2020) is chosen by default. Names mirror datasets hosted on
#' Spectalce (https://oculargeneexpression.org/singlecell)
#'
#' @param label TRUE/FALSE - should the cartoon have text labels for each cell type? True by default.
#'
#' @examples
#' \dontrun{
#' generate_cartoon_data(gene = "RHO", dataset = "godzilla")
#' generate_cartoon_data(gene = "RHO", dataset = "godzilla", label = T, color = "red")
#' }
#'
#' @export
generate_cartoon_data <- function(gene,
                                  dataset = "godzilla",
                                  label = T,
                                  color = "red"){



### 1. Filter the expression matrix for your gene and generate scale

if( ! gene %in% rownames(data_list[["expression_list"]][[dataset]])) stop('requested gene not in this dataset! if this gene is in the GRCh38 notation, try appending -GRCh38 to the end of the gene name')
data_list <- RetinaCartoon::data_list
# data_list the pre-generated list of cartoon outlines and expression matrices. This list
# is pre-generated in data-raw and no adjustment to the parameter should be required. This
# list has three elements (which are themselves lists): 1. a list of expression matrices. 2.
# a list containing cartoon outlines WITHOUT labels. 3. a list containing cartoon outlines WITH
# labels. The character name of the specified dataset is used to refer to the correct expression
# matrix. The user-specified label parameter (see below) specified if the cartoon should include
# labels or not.

expression_df <- data_list[["expression_list"]][[dataset]] %>%
  rownames_to_column("gene_name") %>%
  filter(gene_name == gene) %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("celltype")

colnames(expression_df)[2] <- "expression"

max_expression <- max(expression_df[,2])
max_expression <- max(max_expression, 1)

colfunc <- colorRampPalette(c("grey", color))
my_colors <- colfunc(101)

expression_df <- expression_df %>%
  mutate(percent = round(expression * 100 / max_expression)) %>%
  mutate(color = my_colors[as.integer(percent + 1)])

### 2. Chose the correct xml_object
if(label == TRUE) {
  xml_object <- data_list[["object_labeled_list"]][[dataset]]
} else {
  xml_object <- data_list[["object_unlabeled_list"]][[dataset]]
}


### 3. Generate the indices of each celltype in the current object according to
###    the color scheme below
color_celltype_tibble <- tribble(~hex_color, ~celltype,
                                 "#FFF200", "astrocyte",
                                 "#ED1C24", "muller",
                                 "#00ADEF", "microglia",
                                 "#BF1B2C", "amacrine",
                                 "#BFB290", "horizontal",
                                 "#4F4C4D", "rod",
                                 "#F68712", "bipolar",
                                 "#00A650", "rgc",
                                 "#BDBCBC", "cone",
                                 "#2E3192", "retina-endothelial",
                                 "#EC008C", "retina-pericyte",
                                 "#4AC2C0", "rpe",
                                 "#00ADEF", "t-cell",
                                 "#9C1661", "b-cell",
                                 "#782F9A", "choriocapillaris",
                                 "#E92822", "artery",
                                 "#22275D", "vein",
                                 "#D3AADA", "rpechor-pericyte",
                                 "#F9B0B0", "smc",
                                 "#80CC28", "fibroblast",
                                 "#6F4324", "macrophage",
                                 "#BBE49B", "melanocyte",
                                 "#ED1184", "mast",
                                 "#1D6434", "schwann-nm",
                                 "#F9A72B", "schwann-mye",
                                 "#FCD9E1", "choroid-endothelial"
)

celltype_indexes <- vector(mode = "list", length = nrow(color_celltype_tibble))
names(celltype_indexes) <- color_celltype_tibble[["celltype"]]

all_colors <- c() # grab the color of each path in the object
for(i in 1:length(xml_object@paths)){
  #xml_object@paths[369]$text
  if(! is.null(xml_object@paths[i]$path)){
    color <- xml_object@paths[i]$path@rgb
  } else {
    color <- "BLANK"
  }
  all_colors <- c(all_colors, color)

}


## This for loop returns a list
## Each element in the list corresponds to a different celltype
## Each element in the list is a VECTOR, which corresponds to the path index(es) corresponding to each celltype

for(i in 1:length(celltype_indexes)){
  celltype_index <- which(color_celltype_tibble[["celltype"]] == names(celltype_indexes)[i])
  my_vector <- which(all_colors == color_celltype_tibble[["hex_color"]][celltype_index])
  celltype_indexes[[i]] <- my_vector
}


### 4. Re-color the object based on the color scale from 1 and the indices from 3
for(i in 1:nrow(expression_df)){
  celltype <- expression_df[["celltype"]][i]
  indexes <- unlist(celltype_indexes[[celltype]])
  color <- expression_df[["color"]][i]
  if(length(indexes) > 0){
    for(p in 1:length(indexes)){
      xml_object@paths[indexes[p]]$path@rgb  <- color
    }
  }

}

expression_sequenced <- data.frame(
  expression = seq(from = 0,
                   to = max_expression,
                   by = max_expression/100),
  my_colors)

heatmap_legend <- ggplot(expression_sequenced) +
  geom_tile(aes(x = 1, y = expression, fill = expression)) +
  scale_fill_gradient(low = my_colors[1],
                      high = my_colors[101]) +
  scale_x_continuous(limits=c(0,2),breaks=1, labels = "expression")

leg <- ggpubr::get_legend(heatmap_legend)


return(list(xml_object, leg))

}


#' Plot the generated retina cartoon
#'
#' \code{plot_cartoon} returns a list with two elements. The first element is a
#' heatmap-cartoon genenerated by the grImport package. The second is a legend that related
#' expression to color. These elements can then be combined by the plot_cartoon function
#' to display a pretty cartoon plot
#'
#' @param cartoon_data Cartoon data generated from the generate_cartoon_data
#'
#' @examples
#' \dontrun{
#' generate_cartoon_data(gene = "RHO", dataset = "godzilla") %>% plot_cartoon()
#' }
#'
#' @export
plot_cartoon <- function(cartoon_data){
  cowplot::plot_grid(
    grImport::pictureGrob(cartoon_data[[1]]),
    cartoon_data[[2]],
    ncol = 2,
    rel_widths = c(0.8, 0.2)
  )
}
