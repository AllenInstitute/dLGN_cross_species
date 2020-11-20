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


QC_plot <- function(anno=anno, groupvar="cluster", plot.var="Genes.Detected",stats,statsvar=median) {

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
                  label = stats[[statsvar]]),
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


get_core_transition <- function(norm.dat, cl, select.markers, n.bin=5, n.iter=100, mc.cores=10)
  {
    cl.cv <- parallel::mclapply(1:n.iter, function(i){
      print(i)
      tmp=test_cv_cor(norm.dat, cl, select.markers, n.bin=n.bin)
    }, mc.cores=mc.cores)
    
    cell.cl.cor.map = do.call("rbind",sapply(cl.cv, function(x){
      df = data.frame(cell=names(x),cl=x)
    },simplify=F))
    cell.cl.cor.map = table(cell.cl.cor.map[,1],cell.cl.cor.map[,2])
    cell.cl.cor.map = cell.cl.cor.map / rowSums(cell.cl.cor.map)

    cell.cl.map.df = data.frame(org.cl= as.character(cl[row.names(cell.cl.cor.map)]),best.score=rowMaxs(cell.cl.cor.map), best.cl = colnames(cell.cl.cor.map)[apply(cell.cl.cor.map, 1, which.max)], stringsAsFactors=FALSE)
    row.names(cell.cl.map.df) = row.names(cell.cl.cor.map)
    tmp=get_pair_matrix_coor(cell.cl.cor.map, row.names(cell.cl.map.df), as.character(cell.cl.map.df$best.cl))
    tmp1 = cell.cl.cor.map
    tmp1[tmp]= 0
    cell.cl.map.df$second.score = rowMaxs(tmp1)
    cell.cl.map.df$second.cl =colnames(tmp1)[apply(tmp1,1, which.max)]
    cell.cl.map.df$second.cl[cell.cl.map.df$second.score ==0] = NA
    
    cell.cl.map.df$transition.cl = NA
    tmp = with(cell.cl.map.df, org.cl!=best.cl | best.score < 0.9)
    cell.cl.map.df[tmp,"transition.cl"] = as.character(cell.cl.map.df[tmp,"best.cl"])
    tmp = with(cell.cl.map.df, which(org.cl==transition.cl))
    cell.cl.map.df$transition.cl[tmp] = as.character(cell.cl.map.df[tmp,"second.cl"])
    
    cl.med <- do.call("cbind",tapply(names(cl), cl, function(x){
      rowMedians(as.matrix(norm.dat[select.markers,x]))
    }))
    row.names(cl.med) = select.markers
    
    cell.cl.cor=cor(as.matrix(norm.dat[select.markers, row.names(cell.cl.map.df)]), cl.med[select.markers,])
    cell.cl.map.df$cor = with(cell.cl.map.df, get_pair_matrix(cell.cl.cor, row.names(cell.cl.map.df),as.character(org.cl)))
    cell.cl.map.df$core = is.na(cell.cl.map.df$transition.cl)
    return(cell.cl.map.df)
  }

test_cv_cor <- function(norm.dat, cl, markers, n.bin=5,g.perc=1){
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  }))
  names(bins) = gsub(".*\\.", "", names(bins))
  bins= bins[names(cl)]
  pred.cl = setNames(rep(NA, length(cl)), names(cl))
  for(i in 1:n.bin){
    #print(i)
    train.cells = names(cl)[bins!=i]
    test.cells =names(cl)[bins==i]
    select.markers=sample(markers, round(length(markers)*g.perc))
    map.result <- map_by_cor(norm.dat[select.markers,], cl[train.cells], norm.dat[select.markers, test.cells])$pred.df
    pred.cl[test.cells] = as.character(map.result[test.cells, "pred.cl"])
  }
  return(pred.cl)
}


getConfusionMatrix <- function(orgCluster,
                               foundCluster,
                               proportions = TRUE) {
    orgCluster <- as.character(orgCluster)
    foundCluster <- as.character(foundCluster)
    #lev <- sort(unique(c(orgCluster, foundCluster)))
    #orgCluster <- factor(orgCluster, levels = lev)
    #foundCluster <- factor(orgCluster, levels = lev)
    confusion <- table(foundCluster, orgCluster)
    if (proportions) {
        cs <- colSums(confusion)
        for (i in 1:dim(confusion)[1])
            confusion[i, ] <- confusion[i, ] / pmax(cs, 1e-08)
    }
    confusion
}
