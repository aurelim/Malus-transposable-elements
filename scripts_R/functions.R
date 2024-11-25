window_size = 1e6
max_end = 6e7


# Une fonction pour filtrer les TE avec des coordonnées chevauchantes, en ne conservant que les plus longs 
filter_overlapping_rows_recursive <- function(chromosome_data) {
  # Sort rows by chr and coordinates
  chromosome_data_filtered <- chromosome_data %>%
    arrange(genome, seqid, start) %>%
    mutate(length = end - start + 1) %>%
    
    # Group overlapping rows with the same index
    mutate(cum = head(c(0, cumsum((end < lead(start)) | (end > lead(start) & start > lead(end)))), -1)) %>%
    group_by(cum) %>%
    # Keep only the max length row
    filter(length == max(length)) %>% 
    ungroup() %>%
    select(-cum, -length) 
  
  print(chromosome_data_filtered)
  
  # Check if rows have been removed
  if (nrow(chromosome_data_filtered) < nrow(chromosome_data)) {
    # If yes, recall the function
    return(filter_overlapping_rows_recursive(chromosome_data_filtered))
  } else {
    # If no, return filtered data
    return(chromosome_data)
  }
}

# Une fonction pour ajouter la classification aux TEs
add_classification <- function(chromosome_data, classification_results) {
  chromosome_data_classified <- chromosome_data %>% 
    # mutate(TE = str_remove_all(pattern = "Target [\"]Motif:", attributes)) %>% 
    # mutate(TE = str_remove_all(pattern = "[\"]", TE)) %>% 
    # separate(TE, into = c("TE", "V1", "V2"), sep = " ") %>% 
    # mutate(TE = str_remove_all(TE, "_NA")) %>% 
    # mutate(TE = str_remove_all(TE, "\\)$")) %>% 
    # select(-c(V1, V2, attributes)) %>% 
    left_join(classification_results, by = "TE") %>% 
    select(c(genome, seqid, start, end, strand, TE, Order)) %>% 
    mutate(Order = ifelse(is.na(Order), "Unknown", Order)) %>% 
    mutate(Order = ifelse(Order == "Helitron|SINE", "confused", Order))
}

# Une fonction pour calculer les pourcentage de nucléotides dans chaque fenetre le long de chaque chromosome compris dans un TE, en différenciant les TE par ordre.  
calculate_percentage_positions <- function(chromosome_data, window_size) {
  num_positions_in_table_by_chromosome_and_order <- list()
  
  # Parcourir chaque chromosome unique dans les données
  for (chr in sprintf("chr%s", seq(1:17))) {
    chromosome_subset <- subset(chromosome_data, seqid == chr)
    
    # Parcourir chaque ordre de TE 
    for (order in c("LTR", "TIR", "LINE", "SINE", "noCat")) {
      # for (order in unique(chromosome_data$Order)) {
      # Sous-ensemble de données pour le chromosome actuel et l'ordre de TE actuel
      chromosome_subset_order <- subset(chromosome_subset, Order == order)
      
      # if(nrow(chromosome_subset_order > 0)) {
      
      # Calcul du nombre de positions dans chaque fenêtre par chromosome et ordre
      motif = paste0(chr, "_",order)
      num_positions_in_table_by_chromosome_and_order[[motif]] <- sapply(1:ceiling(max_end / window_size), function(i) {
        window_start <- (i - 1) * window_size
        window_end <- i * window_size
        
        # Sélection des éléments qui chevauchent la fenêtre
        overlapping_elements <- chromosome_subset_order$end >= window_start & chromosome_subset_order$start <= window_end
        
        # Calcul de la longueur pour chaque élément
        element_lengths <- pmin(chromosome_subset_order$end[overlapping_elements], window_end) - pmax(chromosome_subset_order$start[overlapping_elements], window_start) + 1
        
        # Somme des longueurs pour les éléments inclus dans la fenêtre
        sum(element_lengths)
      })
    } 
    
    
    # Délisté les résultats et les fusionner en un seul data frame
    final_result <- do.call(rbind, lapply(names(num_positions_in_table_by_chromosome_and_order), function(motif) {
      data.frame(
        seqid = rep(motif, length(num_positions_in_table_by_chromosome_and_order[[motif]])),
        start_windows = seq(0, by = window_size, length.out = length(num_positions_in_table_by_chromosome_and_order[[motif]])),
        end_windows = seq(window_size, by = window_size, length.out = length(num_positions_in_table_by_chromosome_and_order[[motif]])),
        percentage = num_positions_in_table_by_chromosome_and_order[[motif]] / window_size * 100
      ) %>% 
        separate(seqid, into = c("seqid", "Order"))
      
    }))
    
  }    
  # }
  
  return(final_result)
}

plot_TE_distribution <- function(genome, data_percents, Haut_title) {
  
  df = data_percents %>% 
    dplyr::rename(chr = seqid, start = start_windows, end = end_windows) %>% 
    dplyr::select(chr, start, end, Order, percentage) %>% filter(Order %in% c("noCat", "LTR", "LINE", "SINE", "TIR"))
  
  nObsType = length(unique(df$Order))
  
  df$chr <- factor(df$chr, levels=sprintf("chr%s", seq(1:17)))
  df <- df %>% arrange(chr, start)
  n_windows = ceiling(6e7/1e6)
  
  N = n_windows*17 
  
  df = df %>% mutate(id = rep(seq(1,N), each = nObsType))
  
  
  # Get the name and the y position of each label
  label_data <- df %>% group_by(id, chr) %>% summarize(tot=sum(percentage))
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data <- df %>% 
    group_by(chr) %>% 
    summarize(start=min(id), end=max(id)) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  df$Order <- factor(df$Order, levels = c("noCat", "LTR", "LINE", "SINE", "TIR"))
  
  # Make the plot
  p = df %>% ggplot() +      
    
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=percentage, fill=Order), stat="identity", alpha=1)+
    scale_fill_manual(values = c("LTR" = "orangered4",
                                 "LINE" = "orangered2",
                                 "SINE" = "sienna1",
                                 "TIR" = "blue",
                                 "noCat" = "grey"))+
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    
    ylim(-150,max(label_data$tot, na.rm=T)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(c(-1,0,-1,-1), "cm") 
    ) +
    coord_polar() +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -18, label=chr), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
    ggtitle(label = genome)+
    theme(plot.title = element_text(hjust = 0.5, vjust = Haut_title))
  
  p  
  
  ggsave(plot = p, filename = paste0(genome, "_circularPlot.png"), dpi = 600)
  return(p)
}
