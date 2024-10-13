# Loading of the package used to edit the interface theme
library(fresh)

# Creation of your own Dashboard theme (fresh package).
mytheme <- create_theme(
  
  adminlte_color(
    
    light_blue = "#5D263E"
    
  ),
  
  adminlte_sidebar(
    
    width = "400px",
    dark_bg = "#D8DEE9",
    dark_hover_bg = "#81A1C1",
    dark_color = "#5D263E",
    dark_submenu_color = "#5D263E"
    
  ),
  
  adminlte_global(
    
    content_bg = "#FFF",
    box_bg = "#D8DEE9", 
    info_box_bg = "#D8DEE9"
    
  ), 
  # output_file = "themes/theme_cerulean.css"
  output_file = "themes/theme_orange.css"
)
