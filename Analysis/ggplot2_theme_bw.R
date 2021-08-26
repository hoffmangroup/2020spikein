##ggplot2 theme 
theme_set(theme_bw() +
          theme(axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 16),
      axis.text.y = element_text(vjust = 0.5, size = 16),
      axis.title.x = element_text(vjust = 0.5, size = 20),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      strip.text.x = element_text(size = 16),
      plot.title =
        element_text(hjust = 0.5, color = "black", size = 16, face = "bold")))

