# Load required libraries
library(ggplot2)
library(dplyr)

# Create output directory if needed
if (!dir.exists("output/model_diagrams")) {
  dir.create("output/model_diagrams", recursive = TRUE)
}

# Set random seed for reproducibility
set.seed(42)

# 1. Contemporaneous Model Diagram
# All deposits share a single occupation window
contemp_data <- data.frame(
  deposit_num = c(rep(1, 3), rep(2, 3), rep(3, 2)),
  deposit_label = c(rep("Deposit A", 3), rep("Deposit B", 3), rep("Deposit C", 2)),
  age = c(rnorm(3, 500, 20), rnorm(3, 500, 20), rnorm(2, 500, 20))
)

p1 <- ggplot(contemp_data, aes(x = age, y = deposit_num)) +
  annotate("rect", xmin = 450, xmax = 550, ymin = 0.5, ymax = 3.5,
           fill = "steelblue", alpha = 0.2) +
  geom_point(size = 3, color = "black") +
  annotate("segment", x = 450, xend = 450, y = 0.5, yend = 3.5,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 550, xend = 550, y = 0.5, yend = 3.5,
           color = "red", linewidth = 1.5) +
  annotate("text", x = 500, y = 3.8, label = "Shared occupation window",
           fontface = "italic", size = 4) +
  scale_y_continuous(breaks = c(1, 2, 3),
                     labels = c("Deposit A", "Deposit B", "Deposit C"),
                     limits = c(0.5, 4)) +
  labs(title = "Contemporaneous Model",
       subtitle = "All deposits share a single occupation period",
       x = "Calibrated Age (BP)",
       y = "Archaeological Deposit") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave("output/model_diagrams/model_contemporaneous.png", p1,
       width = 10, height = 6, dpi = 300)

# 2. Sequential Model Diagram
# Deposits are temporally ordered with no overlap
seq_data <- data.frame(
  deposit_num = c(rep(1, 3), rep(2, 3), rep(3, 3)),
  deposit_label = c(rep("Deposit A", 3), rep("Deposit B", 3), rep("Deposit C", 3)),
  age = c(rnorm(3, 400, 15), rnorm(3, 550, 15), rnorm(3, 700, 15))
)

p2 <- ggplot(seq_data, aes(x = age, y = deposit_num)) +
  # Deposit A window
  annotate("rect", xmin = 370, xmax = 430, ymin = 0.7, ymax = 1.3,
           fill = "steelblue", alpha = 0.2) +
  annotate("segment", x = 370, xend = 370, y = 0.7, yend = 1.3,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 430, xend = 430, y = 0.7, yend = 1.3,
           color = "red", linewidth = 1.5) +
  # Deposit B window
  annotate("rect", xmin = 520, xmax = 580, ymin = 1.7, ymax = 2.3,
           fill = "steelblue", alpha = 0.2) +
  annotate("segment", x = 520, xend = 520, y = 1.7, yend = 2.3,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 580, xend = 580, y = 1.7, yend = 2.3,
           color = "red", linewidth = 1.5) +
  # Deposit C window
  annotate("rect", xmin = 670, xmax = 730, ymin = 2.7, ymax = 3.3,
           fill = "steelblue", alpha = 0.2) +
  annotate("segment", x = 670, xend = 670, y = 2.7, yend = 3.3,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 730, xend = 730, y = 2.7, yend = 3.3,
           color = "red", linewidth = 1.5) +
  # Add points
  geom_point(size = 3, color = "black") +
  annotate("text", x = 550, y = 3.6,
           label = "No temporal overlap - ordered succession",
           fontface = "italic", size = 4) +
  scale_y_continuous(breaks = c(1, 2, 3),
                     labels = c("Deposit A", "Deposit B", "Deposit C"),
                     limits = c(0.5, 4)) +
  labs(title = "Sequential Model",
       subtitle = "Deposits are temporally ordered with constraint: end[k] < start[k+1]",
       x = "Calibrated Age (BP)",
       y = "Archaeological Deposit") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave("output/model_diagrams/model_sequential.png", p2,
       width = 10, height = 6, dpi = 300)

# 3. Partial Overlap Model Diagram
# Independent boundaries allowing any overlap pattern
partial_data <- data.frame(
  deposit_num = c(rep(1, 3), rep(2, 3), rep(3, 3)),
  deposit_label = c(rep("Deposit A", 3), rep("Deposit B", 3), rep("Deposit C", 3)),
  age = c(rnorm(3, 500, 20), rnorm(3, 550, 20), rnorm(3, 650, 20))
)

p3 <- ggplot(partial_data, aes(x = age, y = deposit_num)) +
  # Deposit A window
  annotate("rect", xmin = 460, xmax = 540, ymin = 0.7, ymax = 1.3,
           fill = "steelblue", alpha = 0.2) +
  annotate("segment", x = 460, xend = 460, y = 0.7, yend = 1.3,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 540, xend = 540, y = 0.7, yend = 1.3,
           color = "red", linewidth = 1.5) +
  # Deposit B window (overlaps with A)
  annotate("rect", xmin = 510, xmax = 590, ymin = 1.7, ymax = 2.3,
           fill = "steelblue", alpha = 0.2) +
  annotate("segment", x = 510, xend = 510, y = 1.7, yend = 2.3,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 590, xend = 590, y = 1.7, yend = 2.3,
           color = "red", linewidth = 1.5) +
  # Deposit C window (overlaps with B, not A)
  annotate("rect", xmin = 600, xmax = 700, ymin = 2.7, ymax = 3.3,
           fill = "steelblue", alpha = 0.2) +
  annotate("segment", x = 600, xend = 600, y = 2.7, yend = 3.3,
           color = "steelblue", linewidth = 1.5) +
  annotate("segment", x = 700, xend = 700, y = 2.7, yend = 3.3,
           color = "red", linewidth = 1.5) +
  # Add points
  geom_point(size = 3, color = "black") +
  annotate("text", x = 580, y = 3.6,
           label = "Independent boundaries - flexible overlap patterns",
           fontface = "italic", size = 4) +
  scale_y_continuous(breaks = c(1, 2, 3),
                     labels = c("Deposit A", "Deposit B", "Deposit C"),
                     limits = c(0.5, 4)) +
  labs(title = "Partial Overlap Model (Most Flexible)",
       subtitle = "Each deposit has independent start and end dates",
       x = "Calibrated Age (BP)",
       y = "Archaeological Deposit") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave("output/model_diagrams/model_partial_overlap.png", p3,
       width = 10, height = 6, dpi = 300)

cat("\nModel diagrams created successfully!\n")
cat("Output saved to: output/model_diagrams/\n")
cat("  - model_contemporaneous.png\n")
cat("  - model_sequential.png\n")
cat("  - model_partial_overlap.png\n")
