
# Load necessary library
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(VennDiagram)
library(dplyr)
library(VennDiagram)
library(grid)

# Assuming your dataframe is named 'df'
# Count the number of metabolites in each Super_pathway
df<-read.csv("~/DATA overview/cheminfor.csv")

# Assuming your dataframe is named 'df'
# Count the number of metabolites in each Super_pathway
pathway_counts <- df %>%
  group_by(Super_pathway) %>%
  summarise(Count = n())

# Create a horizontal bar plot with customizations
Figure1b_1 <-ggplot(pathway_counts, aes(x = reorder(Super_pathway, Count), y = Count, fill = Super_pathway)) +
  geom_bar(stat = "identity", width = 0.9) +  # Adjust bar width for closer bars
  coord_flip() +  # Flip to make it horizontal
  scale_fill_brewer(palette = "Paired") +  # Use the "Paired" color palette
  geom_text(aes(label = Count), hjust = -0.1, size = 4) +  # Add count labels at the end of bars
  labs(title = "",
       x = "",
       y = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position = "none",  # Remove legend if you prefer
    axis.text.x = element_text(size = 12, color = "black"),  # Customize x-axis text (no bold)
    axis.text.y = element_text(size = 12, color = "black"),  # Customize y-axis text (no bold)
    axis.title.x = element_text(size = 14, color = "black"),  # Customize x-axis label (no bold)
    axis.title.y = element_text(size = 14, color = "black")   # Customize y-axis label (no bold)
  )
Figure1b_1



# Create a data frame with the time points and sample sizes
timepoint_data <- data.frame(
  TimePoint = factor(c("Week 24", "Week 1", "6 Months", "18 Months", "6 Years", "10 Years"),
                     levels = c("Week 24", "Week 1", "6 Months", "18 Months", "6 Years", "10 Years")),  # Factor to preserve order
  Samples = c(581, 579, 520, 535, 477, 520)
)

# Create a bar plot for the number of samples at each time point with enhanced aesthetics
Figure1b_2 <- ggplot(timepoint_data, aes(x = TimePoint, y = Samples)) +
  geom_bar(stat = "identity", fill = "#E9D1D9", color = "black", width = 0.4) +  # Set bar color to grey with black border
  geom_text(aes(label = Samples), vjust = 1.2, size = 5, color = "black") +  # Add sample numbers above the bars
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text.x = element_text(size = 14, color = "black"),  # Customize x-axis text (bold)
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    legend.position = "none",  # Remove legend
    panel.border = element_blank(),  # Remove panel border
    axis.line.x = element_line(color = "black", size = 1),  # Add bold black border to x-axis only
    plot.margin = margin(10, 10, 10, 10)  # Add some padding around the plot for aesthetics
  ) +
  ylim(0, max(timepoint_data$Samples) * 1.1)
Figure1b_2

# Combine the plots
combined_plot <- plot2/ plot1+ plot_layout(heights = c(1, 1.5), guides = "collect") & theme(plot.margin = margin(0, 0, 0, 0))

# Display the combined plot #Figure 1b
combined_plot 


#Venn
load("~/Data/covariat_all_15.Rdata")
load("~/Data/diag_all_name.Rdata")

all_annotation<-read.xlsx("~/Data/Chemical_Annotation_Common.xlsx")
# Add "X" to the start of each cell in the specific column
all_annotation$CHEM_ID <- paste("X", all_annotation$CHEM_ID, sep = "")

longti_data<-read.xlsx("~/Data/all_anchor_norm_metabolome_data.xlsx")
names(longti_data)[4:791] <- paste("X", names(longti_data)[4:791], sep = "")
longti_data<-na.omit(longti_data)
longti_data_1<-longti_data%>%filter(TP%in% c("Mor_w24","Mor_w1","6mth", "18mth", "6yr","10yr"))
#remove unknown
longti_data_2 <- longti_data_1[, !grepl("^X9999", names(longti_data_1))]

filtered_df_diag0<-left_join(diag_all_name[,c(1,3,5,8)],longti_data_2, by="ABCNO")
filtered_df_diag000<-filtered_df_diag0
filtered_df_diag000$newID<- paste(filtered_df_diag000$ABCNO, filtered_df_diag000$TP, sep = "_")

filtered_df_diag1<- filtered_df_diag000[,c(691,1:690)]
filtered_df_diag2<- filtered_df_diag1[!duplicated(filtered_df_diag1$newID), ]

# Assuming your dataframe is called df and has columns: ID, ADHD, NEU_DEV, ASD, TP (for time points)
# Example filtering by time point "Mor_w24"
time_points <- c("Mor_w24", "Mor_w1", "6mth", "18mth", "6yr", "10yr")

# Loop through each time point
for(tp in time_points) {
  # Filter the data for the current time point
  df_filtered <- filtered_df_diag2 %>% filter(TP == tp)
  
  # Create a list of sets for the Venn diagram
  venn_data <- list(
    ADHD = df_filtered$ABCNO[df_filtered$adhd == 1],
    NEU_DEV = df_filtered$ABCNO[df_filtered$neu_dev == 1],
    ASD = df_filtered$ABCNO[df_filtered$asd == 1]
  )
  
  # Customize and plot the Venn diagram
  venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c("ADHD", "NEU_DEV", "ASD"),
    filename = NULL,  # NULL ensures it doesn't write to a file
    output = TRUE,
    fill = c("#E5E5E5", "#D4C0FF", "#C0FFC0"),  # Set custom fill colors for each circle
    alpha = 0.5,  # Set transparency for each circle
    cat.pos = c(-30, 30, 180),  # Adjust position of category labels
    cat.dist = c(0.055, 0.055, 0.055),  # Adjust distance of category labels from the circles
    cat.cex = 2,  # Set the size of category labels
    cex = 1.5,  # Set the size of numbers inside the circles
    fontface = "bold",  # Make the text inside the circles bold
    fontfamily = "sans",  # Use a sans-serif font
    cat.fontface = "bold",  # Make the category labels bold
    cat.fontfamily = "sans",  # Use a sans-serif font for category labels
    lwd = 2,  # Set line width for circle borders
    col = c("black", "black", "black")  # Set border colors for circles
  )
  
  # Display the Venn diagram in the plotting window
  grid.newpage()
  grid.draw(venn.plot)
}
##########four dimension

# Define the time points
time_points <- c("Mor_w24", "Mor_w1", "6mth", "18mth", "6yr", "10yr")

# Loop through each time point
for(tp in time_points) {
  # Filter the data for the current time point
  df_filtered <- filtered_df_diag2 %>% filter(TP == tp)
  
  # Create a list of sets for the Venn diagram
  venn_data <- list(
    ADHD = df_filtered$ABCNO[df_filtered$adhd == 1],
    Any_NDD = df_filtered$ABCNO[df_filtered$neu_dev == 1],
    Autism = df_filtered$ABCNO[df_filtered$asd == 1],
    All = df_filtered$ABCNO  # This includes all samples for the given time point
  )
  
  
  # Customize and plot the 4-set Venn diagram
  venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c("ADHD", "Any_NDD", "Autism", "All"),  # The fourth set represents all samples
    filename = NULL,  # NULL ensures it doesn't write to a file
    output = TRUE,
    
    # Circle properties:
    fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0"),  # Set custom fill colors for circles
    alpha = 0.9,  # Set transparency for circles
    col = "white",  # Set the border color for circles
    lty = 1,  # Line type (solid)
    lwd = 1,  # Set the border width
    
    # Label (intersection values) properties:
    label.col = "black",  # Color for the intersection numbers
    cex = 1.2,  # Font size for intersection numbers
    fontfamily = "serif",  # Font family for intersection numbers
    fontface = "bold",  # Bold the intersection numbers
    
    # Set name (category) properties:
    cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b"),  # Colors for category labels
    cat.cex = 1.5,  # Smaller font size for category labels
    cat.fontfamily = "serif",  # Font family for category labels
    cat.fontface = "bold",  # Make category labels bold
    cat.pos = c(-10, 10, 30, 30),  # Adjust the position of "All" to move it more to the left
    cat.dist = c(0.05, 0.05, 0.05, 0.15),  # Adjust distance of labels from the circles
    cat.default.pos = "outer"  # Ensure labels are outside circles
  )
  
  # Display the Venn diagram in the plotting window
  grid.newpage()
  grid.draw(venn.plot)
  
  # Add the time point label at the top of each plot
  grid.text(tp, y = unit(0.95, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
}

########
# Create a list to store the plots
plots_list <- list()

# Loop through each time point
for(tp in time_points) {
  # Filter the data for the current time point
  df_filtered <- filtered_df_diag2 %>% filter(TP == tp)
  
  # Create a list of sets for the Venn diagram
  venn_data <- list(
    ADHD = df_filtered$ABCNO[df_filtered$adhd == 1],
    Any_NDD = df_filtered$ABCNO[df_filtered$neu_dev == 1],
    Autism = df_filtered$ABCNO[df_filtered$asd == 1],
    All = df_filtered$ABCNO  # This includes all samples for the given time point
  )
  
  # Customize and plot the 4-set Venn diagram
  venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c("ADHD", "Any_NDD", "Autism", "All"),  # The fourth set represents all samples
    filename = NULL,  # NULL ensures it doesn't write to a file
    output = TRUE,
    
    # Circle properties:
    fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0"),  # Set custom fill colors for circles
    alpha = 0.9,  # Set transparency for circles
    col = "white",  # Set the border color for circles
    lty = 1,  # Line type (solid)
    lwd = 1,  # Set the border width
    
    # Label (intersection values) properties:
    label.col = "black",  # Color for the intersection numbers
    cex = 1.2,  # Font size for intersection numbers
    fontfamily = "serif",  # Font family for intersection numbers
    fontface = "bold",  # Bold the intersection numbers
    
    # Set name (category) properties:
    cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b"),  # Colors for category labels
    cat.cex = 0.8,  # Smaller font size for category labels
    cat.fontfamily = "serif",  # Font family for category labels
    cat.fontface = "bold",  # Make category labels bold
    cat.pos = c(-10, 10, 30, 30),  # Adjust the position of "All" to move it more to the left
    cat.dist = c(0.05, 0.05, 0.05, 0.15),  # Adjust distance of labels from the circles
    cat.default.pos = "outer"  # Ensure labels are outside circles
  )
  
  # Create a grob object for the plot
  grid.newpage()
  grid.draw(venn.plot)
  
  # Capture the plot and time point label in one object
  p <- arrangeGrob(
    grobs = list(venn.plot, textGrob(tp, gp = gpar(fontsize = 18, fontface = "bold"))),
    layout_matrix = rbind(c(1), c(2))  # Ensures the label appears above the plot
  )
  
  # Add this plot to the list
  plots_list[[tp]] <- p
}

# Arrange the six Venn diagrams in one row
grid.arrange(grobs = plots_list, nrow = 1)
#############################################
# Assuming your dataframe is called 'df_single'
# Create a list of sets for the Venn diagram based on binary columns
venn_data <- list(
  ADHD = covariat_all_15$ABCNO[covariat_all_15$adhd == 1],
  Any_NDD = covariat_all_15$ABCNO[covariat_all_15$neu_dev == 1],
  Autism = covariat_all_15$ABCNO[covariat_all_15$asd == 1],
  All = covariat_all_15$ABCNO  # This includes all samples
)

# Customize and plot the 4-set Venn diagram
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("ADHD", "Any_NDD", "Autism", "All"),  # The fourth set represents all samples
  filename = NULL,  # NULL ensures it doesn't write to a file
  output = TRUE,
  
  # Circle properties:
  fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0"),  # Set custom fill colors for circles
  alpha = 0.9,  # Set transparency for circles
  col = "white",  # Set the border color for circles
  lty = 1,  # Line type (solid)
  lwd = 1,  # Set the border width
  
  # Label (intersection values) properties:
  label.col = "black",  # Color for the intersection numbers
  cex = 2,  # Font size for intersection numbers
  fontfamily = "serif",  # Font family for intersection numbers
  fontface = "bold",  # Bold the intersection numbers
  
  # Set name (category) properties:
  cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b"),  # Colors for category labels
  cat.cex = 1.5,  # Smaller font size for category labels
  cat.fontfamily = "serif",  # Font family for category labels
  cat.fontface = "bold",  # Make category labels bold
  cat.pos = c(-10, 10, 30, 30),  # Adjust the position of "All" to move it more to the left
  cat.dist = c(0.05, 0.05, 0.05, 0.15),  # Adjust distance of labels from the circles
  cat.default.pos = "outer"  # Ensure labels are outside circles
)

# Display the Venn diagram in the plotting window
grid.newpage()
grid.draw(venn.plot)
###Figure 1b output
