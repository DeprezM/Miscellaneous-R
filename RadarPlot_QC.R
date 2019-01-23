## Marie Deprez -- 23/01/2019
## Radar Plot for scRNA-seq QC metrics

library(tidyverse)
library(scales)
library(ggplot2)
library(ggradar)


# Load csv file with cellRanger QC metrics
sample_stats <- read.csv(file = "RadarPlot_QC_example.csv",
                         dec = ".", sep = ";")


# Format your data
# Add zero minimum values for all QC metrics
data <- sample_stats
data[, c(9:11)] <- data[, c(9:11)] / 100
zero_data <- matrix(0, ncol = 11, nrow = 4)
colnames(zero_data) <- colnames(data)
zero_data[1, 1:2] <- c("Sample_zero_1", "Position_1")
zero_data[2, 1:2] <- c("Sample_zero_2", "Position_2")
zero_data[3, 1:2] <- c("Sample_zero_3", "Position_3")
zero_data[4, 1:2] <- c("Sample_zero_4", "Position_4")

data <- as.data.frame(rbind(data, zero_data))
data[,3:11] <- sapply(data[, 3:11], as.numeric)

# Rescale all values between 0 and 1 (as percentage)
data %>% 
  mutate_at(vars(-Description, -Name, -Sequencing.saturation, 
                 -Fraction.of.Reads.in.Cells, -Reads.Mapped.to.Transcriptome),
            funs(rescale)) -> sample_rescale

sample_rescale <- sample_rescale[grep("zero", sample_rescale$Name, fixed = T, invert = T),]


# plot
ggradar(sample_rescale[sample_rescale$Description == "Position_1", c(1, 3:6, 9:10)], 
        font.radar = "Arial", legend.text.size = 12,
        gridline.mid.colour = "grey",
        axis.label.size = 5, plot.extent.x.sf = 1.75,
        group.line.width = 0.5, axis.label.offset = 1.15,
        group.point.size = 3,grid.label.size = 4) + 
  theme(legend.position = "top")
