#======================================================
# Author: William O. Agyapong
# Purpose: Houses self-defined functions 
# Date created: March 14, 2022
#=======================================================

# This function creates a subgraph for the active subregion including the subregion
# nodes and all other nodes to which there is a connection.
get_subgraph <- function (net, dvwne=F) 
{
  # dvne: delete vertices with no edges
  subregion_df <- subregion_info() # this function is created on the server side
  
  edges_not_needed <- NULL
  
  for (i in 1:length(E(net))) 
  {
    if(!any(subregion_df$ROI %in% unlist(strsplit(as_ids(E(net))[i], "\\|")))) {
      edges_not_needed <- c(edges_not_needed, i)
    }
  }
  
  subg <- delete.edges(net, edges_not_needed)
  if(dvwne == T) {
    subg <- delete.vertices(subg, degree(subg)==0)
  }
  return(subg)
}


#------ This function returns the specified centrality measure from Igraph
# Currently supports only four centrality measures
get_centrality <- function (net, type= c("degree", "closeness", "betweenness", "eigenvector"), 
                            mode=c("sub", "all")) 
{
  type <- match.arg(type)
  mode <- match.arg(mode)
  
  if(type == 'degree') {
    centrality_vals <- degree(net)
    
  } else if (type == "closeness") {
    centrality_vals <- closeness(net)
    
  } else if (type == "betweenness") {
    centrality_vals <- betweenness(net)
    
  } else if (type == "eigenvector") {
    centrality_vals <- eigen_centrality(net)$vector
  }
  
  if(mode == "sub") {
    #---- filter for only the sub-region ----
    subregion_df <- subregion_info() # this function is created on the server side
    return(centrality_vals[subregion_df$ROI] )
  } else {
    return(centrality_vals)
  }
}


#--- custom graph plotting function
# This is a helper function for plotting graphs/networks
pgraph <- function(net, v.size = NULL, dvwne=F, main="", color=NULL, legend=F,
                  toggle.label="none",  edge.col=F, vlabel.size = 0.7, ...)
{
  subregion_df <- subregion_info() # this function is created on the server side
  
  if(dvwne) {
    net <- delete.vertices(net, degree(net)==0)
  }
  
  if(is.null(v.size)) {
    v.size <- 8 # set default vertex size
  } 
  
  if(!is.null(color)) {
    # mark out vertices associated with midbrain
    V(net)$color <- ifelse(V(net)$name %in% subregion_df$ROI, color[1], color[2])
  }
  
  ecol <- ifelse(edge.col, edge_col(net), "gray")
  
  if(toggle.label == "none") {
    vlabel <- NA
  } else if (toggle.label == "short") {
    vlabel <- abbr(V(net)$name)
  } else if (toggle.label == "full") {
    vlabel <- V(net)$name
  } else if (toggle.label == "subgraph") {
    vlabel <- make_vlabel(net)
  }
  
  plot(net,
       vertex.label = vlabel,
       vertex.label.cex = vlabel.size,
       vertex.label.dist=1,
       vertex.size = v.size,
       edge.color = ecol, 
       vertex.label.degree= pi/2,
       # layout = layout_with_fr(net),
       # frame = T
       main = main,
       ...
  )
  if(legend) {
    legend("bottom", legend = c(paste(str_to_title(unique(subregion_df$Subregion)), "subregion"), "Other subregions"),
           pt.cex=1.5, pt.bg = color, pch = 21, cex = 0.7, 
           bty = "n", ncol = 1)
  }
}

# A helper function that helps highlight subregion nodes by making node sizes 
# proportional to the specified centrality measure.
# of the sub-region of interest
make_vsize <- function(subg, centrality = "degree") 
{
  subregion_df <- subregion_info()
  
  centrality_vals <- get_centrality(subg, type = centrality, mode = "all")
  # sub_centrality <- abs(scale(centrality_vals[subregion_df$ROI])[,1])*20
  sub_centrality <- centrality_vals[subregion_df$ROI]
  centrality_new <- vector(length = length(centrality_vals))
  
  # Scale centrality measures for vertices of interest and set others to a fixed constant
  for(i in seq_along(V(subg)$name)) {
    # centrality_new[i] <- ifelse(any(V(subg)$name[i] %in% subregion_df$ROI),
    #                             sub_centrality[names(centrality_vals[i])], 3)
    centrality_new[i] <- ifelse(any(V(subg)$name[i] %in% subregion_df$ROI),
                                centrality_vals[i]/max(sub_centrality)*25, 3)
  }
  
  return(centrality_new)
}


#--- help label only the vertices of interest
make_vlabel <- function(subg) 
{
  subregion_df <- subregion_info()
  
  vnames <- V(subg)$name
  vlabel <- vector(length = length(vnames))
  for(i in seq_along(vnames)) {
    vlabel[i] <- ifelse(any(vnames[i] %in% subregion_df$ROI), abbr(vnames[i]), "")
  }
  return(vlabel)
}

# A special plotting function for subgraphs. Now almost rendered redundant by 
# the pgraph() function above.
psubgraph <- function(subg, centrality = "degree", main="", vlabel.size = 0.9,
                       color=NULL, legend=F, edge.col=F, toggle.label="subgraph") 
{
  
  # vlabel <- ifelse(V(subg)$name %in% subregion_df$ROI, abbr(V(subg)$name), "")
  # vlabel <- make_vlabel(subg)
  vsize <- make_vsize(subg, centrality)
  
  pgraph(subg, 
         vlabel.size = vlabel.size,
         vertex.label.dist=1, 
         vertex.label.color = "darkgreen",
         vertex.label.degree= pi/2,
         v.size = vsize, 
         legend = legend,
         color = color, 
         main =  main,
         toggle.label = toggle.label, 
         edge.curved=F, 
         layout = layout_on_grid(subg))
  
}

#----- This function compares the distributions via a boxplot or density plot
compare_distrib <- function (net1, net2, centrality, plt=c("density", "boxplot", "violin"),
                             annotate.test=F, col_vals = c("gold", "dodgerblue"),
                             method="wilcox", paired=T) 
{ 
  plt <- match.arg(plt) # defaults to density plot
  
  IL_cent <- get_centrality(net1, type = centrality)
  KI_cent <- get_centrality(net2, type = centrality)
  
  df <- data.frame(IL = IL_cent, KI = KI_cent, row.names = NULL) %>%
    pivot_longer(everything(), names_to = "state", values_to = centrality)
  
  test_result <- list(
    if (annotate.test == T) {
      if (method == "wilcox") {
        stat_compare_means(label.sep = " | ", vjust = 1, color = "red",
                           paired = paired)
      } else {
        stat_compare_means(label.sep = " | ", vjust = 1, color = "red",
                           method = "t.test", paired = paired)
      }
    }
    else NULL
  )
  
  if (plt == "boxplot") {
      ggplot(df, aes_string("state", centrality, fill = "state")) +
      geom_boxplot(alpha = 0.3) + xlab("") +
      scale_fill_manual(values = col_vals) +
      theme(legend.position = "none") + test_result
    
  } else if (plt == "density") {
      ggplot(df, aes_string(centrality, fill = "state")) +
      geom_density(alpha = 0.3) + labs(x="", fill = "") +
      scale_fill_manual(values = col_vals) +
      theme_void() + theme(legend.position = "none") + test_result
  } else if (plt == "violin") {
    ggplot(df, aes_string("state", centrality, fill = "state")) +
      geom_violin(alpha = 0.3) +
      geom_boxplot(width = 0.1) + xlab("") +
      scale_fill_manual(values = col_vals) +
      theme(legend.position = "none") + test_result
    
  }
}


#----- statistical test of significance in differences
# Performs comparison tests of differences and returns a table of 
# the test results
test_diff <- function (net1, net2, centrality, test = c("wilcox","ttest"), paired = T, ...) 
{
  test <- match.arg(test)
  IL_cent <- get_centrality(net1, type = centrality)
  KI_cent <- get_centrality(net2, type = centrality)
  
  df <- data.frame(IL = IL_cent, KI = KI_cent, row.names = NULL) %>%
    pivot_longer(everything(), names_to = "state", values_to = "centrality") %>%
    mutate(state = as.factor(state))
  
  if (test == "wilcox") {
    result <- wilcox.test(df$centrality ~ df$state, paired = paired, conf.int =T, ...)
  } else if (test == "ttest") {
    result <- t.test(df$centrality ~ df$state,  paired = paired, conf.int =T, ...)
  }
  
  if (paired) {
    test_descr <- c("wilcox"="Wilcoxon signed rank test", 
                    "ttest"="Paired two sample t-test")
  } else {
    test_descr <- c("wilcox"="Wilcoxon rank sum test", 
                    "ttest"="Two Sample t-test")
  }
  cap_text <- paste(test_descr[test], "for differences in", centrality)
  
  data.frame(statistic = result$statistic, p.value = result$p.value,
             conf.lower = result$conf.int[1], conf.upper = result$conf.int[2]) %>%
    kable(booktabs = T, col.names = c("Test statistic", "P-value", "Conf.lower", "Conf.upper"),
          caption = cap_text, row.names = F) %>%
    kable_paper() %>%
    kable_styling(latex_options = c("HOLD_position"), full_width = F)
}

# Creates acronyms for vertex names
abbr <- function (names.vec, minlength = 3) {
  # names.vec <- subregion_df$ROI
  notation <- vector(length = length(names.vec))
  for(i in seq_along(names.vec)) {
    words <- unlist(str_split(names.vec[i], " "))
    
    if (length(words) == 1 & length(unlist(str_split(words,""))) == minlength) {
      notation[i] <- names.vec[i] # don't abbreviate reasonably short names
    } else if(length(words) < minlength) {
      notation[i] <- abbreviate(str_to_title(names.vec[i]), length(words), method = "both")
    } else  {
      notation[i] <- abbreviate(str_to_title(names.vec[i]), minlength, method = "both")
    }
  }
  return(notation)
  # return(gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1', vnames, perl = TRUE))
}

# color the edges of the graph based on their source node color.
# edge_col <- function(net) {
#   
#   edge.start <- ends(net, es=E(net), names=F)[,1]
#    return(V(net)$color[edge.start])
# }


