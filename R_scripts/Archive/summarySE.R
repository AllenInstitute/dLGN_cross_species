## 
## summarize function
##


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, roundall = F) {
  require(dplyr)
  # This does the summary. For each group return a vector with
  # N, mean, and sd
  
  names(data)[names(data) == measurevar] <- "measurevar"
  
  datac <- data %>%
    select(one_of(groupvars,"measurevar")) %>%
    filter(ifelse(na.rm == T, !is.na(measurevar), T)) %>%
    mutate(measurevar = as.numeric(measurevar)) %>%
    group_by_(c(groupvars)) %>%
    summarise(N = n(),
              median = median(measurevar),
              mean = mean(measurevar),
              max = max(measurevar),
              sd = ifelse(N == 1, 0, sd(measurevar)),
              q25 = as.numeric(quantile(measurevar, 0.25)),
              q75 = as.numeric(quantile(measurevar, 0.75))) %>%
    mutate(se = sd/sqrt(N))
  #%>%
  #  mutate(ci =  se * qt(conf.interval/2 + 0.5, N-1))
  
  
  if(roundall) {
    roundcols <- c("median","mean","max","sd","q25","q75","se","ci")
    datac[roundcols] <- round(datac[roundcols],3)
  }
  
  # datac <- datac %>%
  #   mutate(xpos = 1:n())
  
  return(datac)
}



fix_factors <- function(df) {
  as.data.frame(lapply(df, function(x) {
    if(class(x) == "factor") {
      as.character(x)
    } else {
      x
    }
  }))
}


## 
## qc jitter plot
##


QC_plot <- function(anno=anno, groupvar="cluster", plot.var="Genes.Detected",stats=NULL) {

    groupvar_id <- paste0(groupvar,"_id")
    groupvar_label <- paste0(groupvar,"_label")
    groupvar_color <- paste0(groupvar,"_color")

    max.lim <- 10^ceiling(log10(max(stats$max)))

  plot <- ggplot() +
    # geom_quasirandom from the ggbeeswarm package
    # makes violin-shaped jittered point plots
    geom_quasirandom(data = anno,
                     aes(x = anno[[groupvar_id]],
                         y = anno[[plot.var]]),
                     color = "skyblue",
                     # Need to set position_jitter height = 0 to prevent
                     # jitter on the y-axis, which changes data representation
                     position = position_jitter(width = .3,height = 0),
                     size = 0.1) +
    # Errorbars built using Genes.Detected_stats values
    geom_errorbar(data = stats,
                  aes(x = stats[[groupvar_id]],
                      ymin = q25,
                      ymax = q75),
                  size = 0.2) +
    # Median points from Genes.Detected_stats
    geom_point(data = stats,
               aes(x = stats[[groupvar_id]],
                   y = median),
               color = "red",
               size = 0.5) +
    # Cluster labels as text objects
    geom_text(data = stats,
              aes(x = stats[[groupvar_id]],
                  y = -600,
                  label = stats[[groupvar_label]],
                  color = stats[[groupvar_color]]),
              angle = 90,
              hjust = 1,
              vjust = 0.3,
              size = 3*5/6) +
    # Median values next to cluster labels, since there's space there.
    geom_text(data = stats,
              aes(x = stats[[groupvar_id]],
                  y = -500,
                  label = median),
              angle = 90,
              hjust = 0,
              vjust = 0.3,
              size = 3*5/6) +
    scale_color_identity()+
    # Remove X-axis title
    scale_x_continuous("") +
    theme_bw(6) +
    # Theme tuning
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
    

    return(plot)

}


