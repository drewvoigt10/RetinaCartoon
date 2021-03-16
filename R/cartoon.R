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
#' @param dataset Character string of the dataset name. Options are: "all_retina_rpe_chor", "retina_fov_vs_periph",
#' "retina_AIR_vs_control", "RPE_choroid_unselected", "RPE_choroid_CD31_selected", "retina_fovea_perifovea", and
#' "CD31_choroid_infant_adult". "all_retina_rpe_chor" (a merged object of all retina+rpe+choroid
#' single-cell datasets published by 09/01/2020) is chosen by default. Names mirror datasets hosted on
#' Spectacle (singlecell-eye.org)
#'
#' @param label TRUE/FALSE - should the cartoon have text labels for each cell type? True by default.
#'
#' @examples
#' \dontrun{
#' generate_cartoon_data(gene = "RHO", dataset = "all_retina_rpe_chor")
#' generate_cartoon_data(gene = "RHO", dataset = "all_retina_rpe_chor", label = T, color = "red")
#' }
#'
#' @export
generate_cartoon_data <- function(gene,
                                  dataset = "all_retina_rpe_chor",
                                  label = T,
                                  color = "red"){



### 1. Filter the expression matrix for your gene and generate scale
data_list <- RetinaCartoon::data_list

if( ! gene %in% rownames(data_list[["expression_list"]][[dataset]])) stop('requested gene not in this dataset! if this gene is in the GRCh38 notation, try appending -GRCh38 to the end of the gene name')
# data_list the pre-generated list of cartoon outlines and expression matrices. This list
# is pre-generated in data-raw and no adjustment to the parameter should be required. This
# list has three elements (which are themselves lists): 1. a list of expression matrices. 2.
# a list containing cartoon outlines WITHOUT labels. 3. a list containing cartoon outlines WITH
# labels. The character name of the specified dataset is used to refer to the correct expression
# matrix. The user-specified label parameter (see below) specified if the cartoon should include
# labels or not.

expression_df <- data_list[["expression_list"]][[dataset]] %>%
  tibble::rownames_to_column("gene_name") %>%
  dplyr::filter(gene_name == gene) %>%
  tibble::column_to_rownames("gene_name") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("celltype")

colnames(expression_df)[2] <- "expression"

max_expression <- max(expression_df[,2])
max_expression <- max_expression + 1  # adding 1 to the max expression makes visualizations for
                                      # lowly expressed genes more accurate

colfunc <- grDevices::colorRampPalette(c("grey", color))
my_colors <- colfunc(101)

expression_df <- expression_df %>%
  dplyr::mutate(percent = round(expression * 100 / max_expression)) %>%
  dplyr::mutate(color = my_colors[as.integer(percent + 1)])

### 2. Chose the correct xml_object
if(label == TRUE) {
  xml_object <- data_list[["object_labeled_list"]][[dataset]]
} else {
  xml_object <- data_list[["object_unlabeled_list"]][[dataset]]
}


### 3. Generate the indices of each celltype in the current object according to
###    the color scheme below
color_celltype_tibble <- tibble::tribble(~hex_color, ~celltype,
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

# there are many paths corresponding to different shapes in the xml drawing
# the following for loop searches through all of these paths and extracts the
# hex code color filled in by the path. These colors, in combination with the
# above color_cell type_tibble, will provide us with a way to match each path
# to a cell type

all_path_colors <- rep("", length(xml_object@paths)) # stores the color of each path in the object
for(i in 1:length(xml_object@paths)){

  if(! is.null(xml_object@paths[i]$path)){
    color <- xml_object@paths[i]$path@rgb # extracts the rgb color of the path
  } else {
    color <- "BLANK" # if the path doesnt have a color, we don't care about it!
  }
  all_path_colors[i] <- color
}



# We need to find which paths in our XML object correspond to which shapes.
# First, we create a list that will store the path indexes of each cell type
celltype_path_indexes <- vector(mode = "list", length = nrow(color_celltype_tibble))
names(celltype_path_indexes) <- color_celltype_tibble[["celltype"]]

# Next, we use a for-loop to:
# (1) pull the hexcode color of each cell type of interest
# (2) find which paths in our XML object have this hexcode color

for(i in 1:length(celltype_path_indexes)){

  celltype_color <- color_celltype_tibble %>%
    dplyr::filter(celltype == names(celltype_path_indexes)[i]) %>%
    dplyr::pull(hex_color)

  path_index <- which(all_path_colors == celltype_color)
  celltype_path_indexes[[i]] <- path_index
}


### 4. Re-color the object based on the EXPRESSION color scale
# (re-coloring modifies the original XML object)
for(i in 1:nrow(expression_df)){
  celltype <- expression_df[["celltype"]][i]
  indexes <- unlist(celltype_path_indexes[[celltype]])
  color <- expression_df[["color"]][i]
  if(length(indexes) > 0){
    for(p in 1:length(indexes)){
      xml_object@paths[indexes[p]]$path@rgb  <- color
    }
  }

}

### 5. Creating a legend of expression
expression_sequenced <- data.frame(
  expression = seq(from = 0,
                   to = max_expression,
                   by = max_expression/100),
  my_colors)

heatmap_legend <- ggplot2::ggplot(expression_sequenced) +
  ggplot2::geom_tile(aes(x = 1, y = expression, fill = expression)) +
  ggplot2::scale_fill_gradient(low = my_colors[1],
                      high = my_colors[101]) +
  ggplot2::scale_x_continuous(limits=c(0,2),breaks=1, labels = "expression")

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
#' @param cartoon_data Cartoon data generated from the generate_cartoon_data. Cartoon data is
#' a list containing an XML object and an expression legend.
#'
#' @examples
#' \dontrun{
#' generate_cartoon_data(gene = "RHO", dataset = "all_retina_rpe_chor") %>% plot_cartoon()
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
