rm(list=ls(all=TRUE)) #Clear the memory of variables from previous run. This is not called by knitr, because it's above the first chunk.

# ---- load_sources ------------------------------------------------------------
source("http://statpower.net/Content/312/R%20Stuff/Steiger%20R%20Library%20Functions.txt")
source("https://raw.githubusercontent.com/andkov/psy532/master/scripts/graphs/main_theme.R")

# ---- load_packages -----------------------------------------------------------
# install.packages(c("vegan", "ecodist", "labdsv", "ape", "ade4", "smacof"))
library(ggplot2) #For graphing
library(magrittr) #Pipes
requireNamespace("knitr", quietly=TRUE)
requireNamespace("scales", quietly=TRUE) #For formating values in graphs
# requireNamespace("RColorBrewer", quietly=TRUE)
requireNamespace("dplyr", quietly=TRUE)
# requireNamespace("plyr", quietly=TRUE)
# requireNamespace("reshape2", quietly=TRUE) #For converting wide to long
# requireNamespace("mgcv, quietly=TRUE) #For the Generalized Additive Model that smooths the longitudinal graphs.
requireNamespace("vegan", quietly = TRUE)
requireNamespace("ecodist", quietly = TRUE)
requireNamespace("labdsv", quietly = TRUE)
requireNamespace("ape", quietly = TRUE)
requireNamespace("ade4", quietly = TRUE)
requireNamespace("smacof", quietly = TRUE)
# ---- declare_globals ---------------------------------------------------------
options(show.signif.stars=F) #Turn off the annotations on p-values

path_input <- "./data-phi-free/derived/motor-trend-car-test.rds"

histogram_discrete <- function(
  d_observed,
  variable_name,
  levels_to_exclude   = character(0),
  main_title          = variable_name,
  x_title             = NULL,
  y_title             = "Number of Included Records",
  text_size_percentage= 6,
  bin_width           = 1L) {

  d_observed <- as.data.frame(d_observed) #Hack so dplyr datasets don't mess up things
  if( !base::is.factor(d_observed[, variable_name]) )
    d_observed[, variable_name] <- base::factor(d_observed[, variable_name])

  d_observed$iv <- base::ordered(d_observed[, variable_name], levels=rev(levels(d_observed[, variable_name])))

  ds_count <- plyr::count(d_observed, vars=c("iv"))
  # if( base::length(levels_to_exclude)>0 ) { }
  ds_count <- ds_count[!(ds_count$iv %in% levels_to_exclude), ]

  ds_summary <- plyr::ddply(ds_count, .variables=NULL, transform, count=freq, proportion = freq/sum(freq) )
  ds_summary$percentage <- base::paste0(base::round(ds_summary$proportion*100), "%")

  y_title <- base::paste0(y_title, " (n=", scales::comma(base::sum(ds_summary$freq)), ")")

  g <- ggplot(ds_summary, aes_string(x="iv", y="count", fill="iv", label="percentage")) +
    geom_bar(stat="identity") +
    geom_text(stat="identity", size=text_size_percentage, hjust=.8) +
    scale_y_continuous(labels=scales::comma_format()) +
    labs(title=main_title, x=x_title, y=y_title) +
    coord_flip()

  theme  <- theme_light(base_size=14) +
    theme(legend.position = "none") +
    theme(axis.text.x=element_text(colour="gray40")) +
    theme(axis.title.x=element_text(colour="gray40")) +
    theme(axis.text.y=element_text(size=14)) +
    theme(panel.border = element_rect(colour="gray80")) +
    theme(axis.ticks.length = grid::unit(0, "cm"))

  return( g + theme )
}
histogram_continuous <- function(
  d_observed,
  variable_name,
  bin_width      = NULL,
  main_title     = variable_name,
  x_title        = paste0(variable_name, " (each bin is ", scales::comma(bin_width), " units wide)"),
  y_title        = "Frequency",
  rounded_digits = 0L
  ) {

  d_observed <- as.data.frame(d_observed) #Hack so dplyr datasets don't mess up things
  d_observed <- d_observed[!base::is.na(d_observed[, variable_name]), ]

  ds_mid_points <- base::data.frame(label=c("italic(X)[50]", "bar(italic(X))"), stringsAsFactors=FALSE)
  ds_mid_points$value <- c(stats::median(d_observed[, variable_name]), base::mean(d_observed[, variable_name]))
  ds_mid_points$value_rounded <- base::round(ds_mid_points$value, rounded_digits)

  g <- ggplot(d_observed, aes_string(x=variable_name)) +
    geom_bar(stat="bin", binwidth=bin_width, fill="gray70", color="gray90", position=position_identity()) +
    geom_vline(xintercept=ds_mid_points$value, color="gray30") +
    geom_text(data=ds_mid_points, aes_string(x="value", y=0, label="value_rounded"), color="tomato", hjust=c(1, 0), vjust=.5) +
    scale_x_continuous(labels=scales::comma_format()) +
    scale_y_continuous(labels=scales::comma_format()) +
    labs(title=main_title, x=x_title, y=y_title) +
    theme_light() +
    theme(axis.ticks.length = grid::unit(0, "cm"))

  ds_mid_points$top <- stats::quantile(ggplot2::ggplot_build(g)$panel$ranges[[1]]$y.range, .8)
  g <- g + ggplot2::geom_text(data=ds_mid_points, ggplot2::aes_string(x="value", y="top", label="label"), color="tomato", hjust=c(1, 0), parse=TRUE)
  return( g )
}
tonics <- c("Canada Dry","Fever Tree", "Great Value",  "Schweppes", "Q"  )

baseSize <- 10
main_theme <- theme_bw() +
  theme(axis.text = element_text(colour="gray40")) +
  theme(axis.title = element_text(colour="gray40")) +
  theme(panel.border = element_rect(colour="gray80")) +
  theme(axis.ticks = element_line(colour="gray80"))

# ---- load_data ---------------------------------------------------------------


  
# ---- tweak_data --------------------------------------------------------------


# ---- define_mds_plotting --------------------------------------

mds_scatter <- function(input, limit){
  dsmat <- CompleteSymmetricMatrix(input)
  (rownames(dsmat) <- tonics)
  (colnames(dsmat) <- tonics)
  (ds <- as.dist(dsmat)  )
  (mds <- cmdscale(ds, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE))
  (d <- as.data.frame(mds))
  
  g <- ggplot2::ggplot(d, aes(x=V1, y=V2)) +
    geom_point(shape=21, fill="red", size=4) +
    geom_text(aes(label=tonics), vjust=-1.5, size=8) +
    scale_y_continuous(limits=c(-limit,limit))+
    scale_x_continuous(limits=c(-limit,limit))+
    main_theme
  return(g)
}

# mds_scatter(input, 5)
