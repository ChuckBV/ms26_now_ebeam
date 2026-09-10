theme_fullwidth <- function(){
  theme_bw() +
    theme(
      axis.text.x = element_text(color = "black", size = 9),#, angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black", size = 9),
      axis.title.x = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 12),
      legend.title = element_text(color = "black", size = 9),
      legend.text = element_text(color = "black", size = 8))
}

theme_powerpoint <- function(){
  theme_bw() +
    theme(
      axis.text.x = element_text(color = "black", size = 12),#, angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 10))
}