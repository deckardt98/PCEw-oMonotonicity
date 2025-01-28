### This file contains simulation code for the article entitled "Quadruply robust 
### principal stratification analysis without monotonicity" by Tong et al..
### SACE

library(truncnorm)
library(SuperLearner)
library(randomForest)
library(geex)
library(parallel)
library(data.table)
library(boot)
library(caret)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)


setwd("~/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/Code")
source("EST_FUN11.R")
source("EST_FUN01.R")
source("EST_FUN00.R")
source("EST_FUN10.R")

#Function simulating full data vector
simu_full_data <- function(n, theta){
  #Input:
  #n: sample size
  #theta: odds ratio
  
  Z <- D <- c()
  #simulate covariates
  X1 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X2 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X3 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X4 <- rbinom(n, size = 1, prob = 0.5)
  if (theta>9999){
    X1 <- abs(X1)
    X2 <- abs(X2)
    X3 <- abs(X3)
    probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
    #max(abs(probZ-0.5))
    Z <- rbinom(n, size = 1, prob = probZ)
    #mean(Z)
    
    #simulate G=(D(0),D(1))
    probD1 <- expit(1.2*X4)
    probD0 <- expit(-0.4-0.2*X1-0.2*X2-0.2*X3-0.2*X4)
    prob11 <- probD0
    prob01 <- probD1 - probD0
    prob00 <- 1 - probD1
    prob_matrix <- cbind(prob00, prob01, prob11)
    # Simulate G from the categorical distribution
    G <- apply(prob_matrix, 1, function(p) sample(0:2, size = 1, prob = p))
    #simulate D(1)
    D1 <- as.numeric(I(G!=0))
    #simulate D(0)
    D0 <- as.numeric(I(G==2))
    #compute observed intermediate outcome
    D <- D1*Z+(1-Z)*D0
  } else {
    probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
    #max(abs(probZ-0.5))
    Z <- rbinom(n, size = 1, prob = probZ)
    #mean(Z)
    #simulate (D(0),D(1))
    probD1 <- expit(0.3*X1+0.4*X2+0.3*X3+0.5*X4)
    probD0 <- expit(0.4*X1+0.3*X2+0.4*X3+0.5*X4)
    deltaX <- (1+(theta-1)*(probD1+probD0))^2-4*theta*(theta-1)*probD1*probD0
    prob11 <- (1+(theta-1)*(probD1+probD0)-sqrt(deltaX))/2/(theta-1)
    prob10 <- probD0 - prob11
    prob01 <- probD1 - prob11
    prob00 <- 1-prob11-prob10-prob01
    prob_matrix <- cbind(prob10, prob00, prob01, prob11)
    jointD <- apply(prob_matrix, 1, function(p) sample(1:4, size = 1, prob = p))
    #simulate D(0)
    D0 <- as.numeric(I(jointD==1|jointD==4))
    #simulate D(1)
    D1 <- as.numeric(I(jointD==3|jointD==4))
    #compute observed intermediate outcome
    D <- D1*Z+(1-Z)*D0
    #max(max(abs(0.5-probD1)),max(abs(0.5-probD0)))
    #table(D1,D0)/n
  }
  #simulate potential final outcome
  meanY1 <- -1+D1+X1+3*X2+3*X3+3*X4
  meanY0 <- 3-D0-1.5*X1+2*X2+2*X3-2*X4
  Y1 <- rnorm(n = n, mean = meanY1, sd = 1)
  Y0 <- rnorm(n = n, mean = meanY0, sd = 1)
  Y <- Y1*Z+(1-Z)*Y0
  
  return(list(X1, X2, X3, X4, Z, D1, D0, D, Y1, Y0, Y))
}


#Approximate truth
cal_tr <-  function(n, theta){
  full_data <-  simu_full_data(n, theta=theta)
  X1 <- full_data[[1]]
  X2 <- full_data[[2]]
  X3 <- full_data[[3]]
  X4 <- full_data[[4]]
  Z <- full_data[[5]]
  D1 <- full_data[[6]]
  D0 <- full_data[[7]]
  D <- full_data[[8]]
  Y1 <- full_data[[9]]
  Y0 <- full_data[[10]]
  Y <- full_data[[11]]
  #approximate truth
  truepce11 <- mean(D0*D1*(Y1-Y0))/mean(D0*D1)
  truepce01 <- mean((1-D0)*D1*(Y1-Y0))/mean((1-D0)*D1)
  truepce00 <- mean((1-D0)*(1-D1)*(Y1-Y0))/mean((1-D0)*(1-D1))
  truepce10 <- mean(D0*(1-D1)*(Y1-Y0))/mean(D0*(1-D1))
  return(c(truepce11,truepce01,truepce00,truepce10))
}


sim_bias <-  function(n, B, truepce, theta){
  bias11 <- bias01 <- bias00 <- bias10 <- matrix(NA, nrow = B, ncol = 41)
  aese11 <- aese01 <- aese00 <- aese10 <- matrix(NA, nrow = B, ncol = 41)
  cov11 <- cov01 <- cov00 <- cov10 <- matrix(NA, nrow = B, ncol = 41)
  for (j in 1:B){
    tryCatch({
      full_data <-  simu_full_data(n, theta)
      X1 <- full_data[[1]]
      X2 <- full_data[[2]]
      X3 <- full_data[[3]]
      X4 <- full_data[[4]]
      Z <- full_data[[5]]
      D <- full_data[[8]]
      Y <- full_data[[11]]
      df <- as.data.frame(cbind(X1,X2,X3,X4,Z,D,Y))
      start <- Sys.time()
      result_est11 <- cal_estimates11(df, theta, truepce[1])
      result_est01 <- cal_estimates01(df, theta, truepce[2])
      result_est00 <- cal_estimates00(df, theta, truepce[3])
      result_est10 <- cal_estimates10(df, theta, truepce[4])
      print( Sys.time() - start )
      bias11[j,] <- result_est11[seq(1,82, by =2)]
      aese11[j,] <- result_est11[-seq(1,82, by =2)]
      bias01[j,] <- result_est01[seq(1,82, by =2)]
      aese01[j,] <- result_est01[-seq(1,82, by =2)]
      bias00[j,] <- result_est00[seq(1,82, by =2)]
      aese00[j,] <- result_est00[-seq(1,82, by =2)]
      bias10[j,] <- result_est10[seq(1,82, by =2)]
      aese10[j,] <- result_est10[-seq(1,82, by =2)]
      #print(round(aese[j,],2))
      cov11[j,] <- I(truepce[1]>=(bias11[j,]-aese11[j,]*1.96))*I(truepce[1]<=(bias11[j,]+aese11[j,]*1.96))
      cov01[j,] <- I(truepce[2]>=(bias01[j,]-aese01[j,]*1.96))*I(truepce[2]<=(bias01[j,]+aese01[j,]*1.96))
      cov00[j,] <- I(truepce[3]>=(bias00[j,]-aese00[j,]*1.96))*I(truepce[3]<=(bias00[j,]+aese00[j,]*1.96))
      cov10[j,] <- I(truepce[4]>=(bias10[j,]-aese10[j,]*1.96))*I(truepce[4]<=(bias10[j,]+aese10[j,]*1.96))
    }, 
    warning = function(w){
      message("Warning in iteration ", j, ": ", conditionMessage(w))
      return(NULL)
    },
    error = function(e) {
      message("Error in iteration ", j, ": ", conditionMessage(e))
      return(NULL)})
  }
  
  return(list(bias11,aese11,cov11,
              bias01,aese01,cov01,
              bias00,aese00,cov00,
              bias10,aese10,cov10))
}


set.seed(920784642)
theta <- 999999
truepce <- cal_tr(10000000, theta)
n <- 500
B <- 1000
result <- sim_bias(n,B,truepce, theta)
copy_result <- result
mcsd11 <- round(apply(result[[1]], 2, sd, na.rm = TRUE),3)
mcsd01 <- round(apply(result[[4]], 2, sd, na.rm = TRUE),3)
mcsd00 <- round(apply(result[[7]], 2, sd, na.rm = TRUE),3)
mcsd10 <- round(apply(result[[10]], 2, sd, na.rm = TRUE),3)
#mcsd
mbias11 <- round(colMeans(result[[1]], na.rm = TRUE)-truepce[1],3)
mbias01 <- round(colMeans(result[[4]], na.rm = TRUE)-truepce[2],3)
mbias00 <- round(colMeans(result[[7]], na.rm = TRUE)-truepce[3],3)
mbias10 <- round(colMeans(result[[10]], na.rm = TRUE)-truepce[4],3)
#mbias
aese11 <- round(colMeans(result[[2]], na.rm = TRUE),3)
aese01 <- round(colMeans(result[[5]], na.rm = TRUE),3)
aese00 <- round(colMeans(result[[8]], na.rm = TRUE),3)
aese10 <- round(colMeans(result[[11]], na.rm = TRUE),3)
#aese
cp11 <- round(colMeans(result[[3]], na.rm = TRUE),3)
cp01 <- round(colMeans(result[[6]], na.rm = TRUE),3)
cp00 <- round(colMeans(result[[9]], na.rm = TRUE),3)
cp10 <- round(colMeans(result[[12]], na.rm = TRUE),3)
#cp

out_table11 <- as.data.frame(cbind(mbias11,mcsd11,aese11,cp11))
out_table01 <- as.data.frame(cbind(mbias01,mcsd01,aese01,cp01))
out_table00 <- as.data.frame(cbind(mbias00,mcsd00,aese00,cp00))
out_table10 <- as.data.frame(cbind(mbias10,mcsd10,aese10,cp10))
out_table <- cbind(out_table11,out_table01,out_table00,out_table10)
rownames(out_table) <- c("Weighting",
                         "PAR0.5_1", "PAR1_1", "PAR2_1", "PARinfty_1", 
                         "ML0.5_1", "ML1_1", "ML2_1", "MLinfty_1", 
                         "PAR0.5_2", "PAR1_2", "PAR2_2", "PARinfty_2", 
                         "ML0.5_2", "ML1_2", "ML2_2", "MLinfty_2", 
                         "PAR0.5_3", "PAR1_3", "PAR2_3", "PARinfty_3", 
                         "ML0.5_3", "ML1_3", "ML2_3", "MLinfty_3", 
                         "PAR0.5_4", "PAR1_4", "PAR2_4", "PARinfty_4", 
                         "ML0.5_4", "ML1_4", "ML2_4", "MLinfty_4",
                         "PAR0.5_5", "PAR1_5", "PAR2_5", "PARinfty_5", 
                         "ML0.5_5", "ML1_5", "ML2_5", "MLinfty_5")
out_table
setwd("/Users/deckard/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/results")
saveRDS(copy_result, "table_new_1.rds")
write.csv(out_table, file= "table_new_1.csv")

#result <- readRDS(file="table_new_0.5.rds")

##### Box plot old #########
# Define a function to process data and create plots
create_bias_plot <- function(design_data, true_value, plot_title, whichcol, est_name, est_color) {
  bias_data <- design_data[,whichcol] - true_value
  df_plot <- as.data.frame(bias_data)
  colnames(df_plot) <- est_name
  df_plot$row_id <- 1:nrow(df_plot)
  df_long <- melt(df_plot, id.vars = "row_id", variable.name = "Column", value.name = "Value")
  
  # Create and return the box plot
  ggplot(df_long, aes(x = Column, y = Value, fill = Column)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "", y = plot_title) +
    coord_cartesian(ylim = c(-1.6, 1.6)) +
    scale_fill_manual(values = est_color) +
    geom_hline(yintercept=0, linetype="dashed")+
    theme(
      axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none",
      panel.background = element_blank(),     # Remove background
      panel.grid.major = element_blank(),     # Remove major grid lines
      panel.grid.minor = element_blank(),     # Remove minor grid lines
      axis.line = element_line(color = "black"),  # Keep axis lines
      axis.ticks = element_line(color = "black")  # Keep axis ticks
    )
}

result[[10]][,c(5,9,13,17,21,25,29,33,37,41)] <- NA
## Design matrix 1
# Plot for 11
box11_1 <- create_bias_plot(result[[1]], truepce[1], "(i)", c(2:9), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                       "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray"))
# Plot for 01 
box01_1 <-  create_bias_plot(result[[4]], truepce[2], "(i)", c(2:9), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                        "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                              c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
# Plot for 00
box00_1 <- create_bias_plot(result[[7]], truepce[3], "(i)", c(2:9), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                               "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray")) 
# Plot for 10
box10_1 <- create_bias_plot(result[[10]], truepce[4], "(i)", c(2:9), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                         "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "orange", "gray", "gray", "orange"))  + theme(axis.title.y = element_blank(), axis.text.y = element_blank())

## Design matrix 2
# Plot for 11
box11_2 <- create_bias_plot(result[[1]], truepce[1], "(ii)", c(10:17), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                         "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray"))
# Plot for 01 
box01_2 <-  create_bias_plot(result[[4]], truepce[2], "(ii)", c(10:17), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                          "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                              c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
# Plot for 00
box00_2 <- create_bias_plot(result[[7]], truepce[3], "(ii)", c(10:17), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                                 "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray"))
# Plot for 10
box10_2 <- create_bias_plot(result[[10]], truepce[4], "(ii)", c(10:17), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                          "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "orange")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())

## Design matrix 3
# Plot for 11
box11_3 <- create_bias_plot(result[[1]], truepce[1], "(iii)", c(18:25), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                         "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray"))
# Plot for 01 
box01_3 <-  create_bias_plot(result[[4]], truepce[2], "(iii)", c(18:25), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                          "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                              c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
# Plot for 00
box00_3 <- create_bias_plot(result[[7]], truepce[3], "(iii)", c(18:25), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                                 "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "gray", "orange", "gray"))
# Plot for 10
box10_3 <- create_bias_plot(result[[10]], truepce[4], "(iii)", c(18:25), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                           "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "orange", "gray", "gray", "orange")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())


## Design matrix 4
# Plot for 11
box11_4 <- create_bias_plot(result[[1]], truepce[1], "(iv)", c(26:33), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                         "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "gray", "gray", "gray", "gray", "orange", "gray")) 
# Plot for 01 
box01_4 <-  create_bias_plot(result[[4]], truepce[2], "(iv)", c(26:33), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                          "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "gray", "gray", "gray", "gray", "orange", "gray"))  + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
# Plot for 00
box00_4 <- create_bias_plot(result[[7]], truepce[3], "(iv)", c(26:33), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                                 "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "gray", "gray", "gray", "gray", "orange", "gray")) 
# Plot for 10
box10_4 <- create_bias_plot(result[[10]], truepce[4], "(iv)", c(26:33), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                          "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "gray", "gray", "gray", "orange")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())

## Design matrix 5
# Plot for 11
box11_5 <- create_bias_plot(result[[1]], truepce[1], "(v)", c(34:41), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                         "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "gray", "gray", "gray", "gray", "orange", "gray")) 
# Plot for 01 
box01_5 <-  create_bias_plot(result[[4]], truepce[2], "(v)", c(34:41), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                          "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                             c("gray", "gray", "gray", "gray", "gray", "gray", "orange", "gray"))  + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
# Plot for 00
box00_5 <- create_bias_plot(result[[7]], truepce[3], "(v)", c(34:41), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                                 "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "gray", "gray", "gray", "gray", "orange", "gray")) 
# Plot for 10
box10_5 <- create_bias_plot(result[[10]], truepce[4], "(v)", c(34:41), c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                                                                         "ML-0.5", "ML-1", "ML-2", "ML-INF"),
                            c("gray", "gray", "gray", "gray", "gray", "orange")) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())


# Convert the ggplot objects to grobs
box11_1_grob <- ggplotGrob(box11_1)
box01_1_grob <- ggplotGrob(box01_1)
box11_2_grob <- ggplotGrob(box11_2)
box01_2_grob <- ggplotGrob(box01_2)
box11_3_grob <- ggplotGrob(box11_3)
box01_3_grob <- ggplotGrob(box01_3)
box11_4_grob <- ggplotGrob(box11_4)
box01_4_grob <- ggplotGrob(box01_4)
box11_5_grob <- ggplotGrob(box11_5)
box01_5_grob <- ggplotGrob(box01_5)

# Create text grobs for the column titles
col1_title <- textGrob("Always-takers", gp = gpar(fontsize = 14, fontface = "bold"))
col2_title <- textGrob("Compliers", gp = gpar(fontsize = 14, fontface = "bold"))

# save pdf
pdf("box2_11_01.pdf", width = 10, height = 8)  # Specify the size here
# Arrange the column titles and the plots
grid.arrange(
  # Combine the column titles in a row
  arrangeGrob(col1_title, col2_title, ncol = 2),
  # Combine the plots in the grid
  arrangeGrob(box11_1_grob, box01_1_grob,
              box11_2_grob, box01_2_grob,
              box11_3_grob, box01_3_grob,
              box11_4_grob, box01_4_grob,
              box11_5_grob, box01_5_grob, 
              nrow = 5, ncol = 2),
  # Specify the number of rows (1 row for titles, 1 for the grid)
  nrow = 2, heights = c(0.1, 1)
)
dev.off()  # Close the PDF device


box00_1_grob <- ggplotGrob(box00_1)
box00_2_grob <- ggplotGrob(box00_2)
box00_3_grob <- ggplotGrob(box00_3)
box00_4_grob <- ggplotGrob(box00_4)
box00_5_grob <- ggplotGrob(box00_5)

box10_1_grob <- ggplotGrob(box10_1)
box10_2_grob <- ggplotGrob(box10_2)
box10_3_grob <- ggplotGrob(box10_3)
box10_4_grob <- ggplotGrob(box10_4)
box10_5_grob <- ggplotGrob(box10_5)

col1_title <- textGrob("Never-takers", gp = gpar(fontsize = 14, fontface = "bold"))
col2_title <- textGrob("Defiers", gp = gpar(fontsize = 14, fontface = "bold"))

pdf("box2_00_10.pdf", width = 10, height = 8)  # Specify the size here
# Arrange the column titles and the plots
grid.arrange(
  # Combine the column titles in a row
  arrangeGrob(col1_title, col2_title, ncol = 2),
  # Combine the plots in the grid
  arrangeGrob(box00_1_grob, box10_1_grob,
              box00_2_grob, box10_2_grob,
              box00_3_grob, box10_3_grob,
              box00_4_grob, box10_4_grob,
              box00_5_grob, box10_5_grob, 
              nrow = 5, ncol = 2),
  # Specify the number of rows (1 row for titles, 1 for the grid)
  nrow = 2, heights = c(0.1, 1)
)
dev.off()  # Close the PDF device

############################


###### New plots #####
# Create a combined dataframe for all plots
setwd("/Users/deckard/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/results")
result <- readRDS(file="table_new_infty.rds")
combined_df <- data.frame()

## 11 01
# Add condition labels
conditions <- c("11", "01")  # Left and right panels
conditions <- factor(conditions, levels = c("11", "01"))
designs <- c("(i)", "(ii)", "(iii)", "(iv)", "(v)")  # Row labels

for (d in 1:5) {
  for (c in 1:2) {
    # Select appropriate result based on condition
    result_idx <- if(c == 1) 1 else 4  # Always-takers = 1, Compliers = 4
    
    # Calculate column indices for the design
    col_start <- 2 + (d-1)*8
    col_end <- col_start + 7
    
    # Extract and process data
    bias_data <- result[[result_idx]][,col_start:col_end] - truepce[if(c == 1) 1 else 2]
    df_temp <- as.data.frame(bias_data)
    colnames(df_temp) <- c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                           "ML-0.5", "ML-1", "ML-2", "ML-INF")
    
    # Add metadata columns
    df_temp$design <- designs[d]
    df_temp$condition <- conditions[c]
    df_temp$row_id <- 1:nrow(df_temp)
    
    # Determine which estimators should be orange
    if (conditions[c] == "Always-takers") {
      if (designs[d] %in% c("(i)", "(ii)", "(iii)")) {
        highlight_est <- c("PAR-INF", "ML-INF")
      } else {
        highlight_est <- c("ML-INF")
      }
    } else {  # Compliers
      if (designs[d] %in% c("(i)", "(ii)", "(iii)")) {
        highlight_est <- c("PAR-INF", "ML-INF")
      } else {
        highlight_est <- c("ML-INF")
      }
    }
    
    # Reshape to long format
    df_long <- reshape2::melt(df_temp, 
                              id.vars = c("row_id", "design", "condition"),
                              variable.name = "estimator",
                              value.name = "bias")
    
    # Add highlight column
    df_long$highlight <- df_long$estimator %in% highlight_est
    
    combined_df <- rbind(combined_df, df_long)
  }
}

### save pdf
pdf("boxINF_11_01.pdf", width = 10, height = 8)  # Specify the size here
# Create the plot
ggplot(combined_df, aes(x = estimator, y = bias, fill = highlight)) +
  geom_boxplot() +
  facet_grid(design ~ condition) +
  coord_cartesian(ylim = c(-1.6, 1.6)) +
  scale_fill_manual(values = c("gray", "orange")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "lightblue", color = "black"),  # Set background color
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    axis.title = element_text(size = 14)
  ) +
  labs(x = "", y = "Bias")
dev.off()  # Close the PDF device

#####

## 00 10
# Create a combined dataframe for all plots
setwd("/Users/deckard/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/results")
result <- readRDS(file="table_new_infty.rds")
result[[10]][,c(5,9,13,17,21,25,29,33,37,41)] <- NA
combined_df <- data.frame()
# Add condition labels
conditions <- c("00", "10")  # Left and right panels
conditions <- factor(conditions, levels = c("00", "10"))
designs <- c("(i)", "(ii)", "(iii)", "(iv)", "(v)")  # Row labels

for (d in 1:5) {
  for (c in 1:2) {
    # Select appropriate result based on condition
    result_idx <- if(c == 1) 7 else 10  # 00 = 7, 10 = 10
    
    # Calculate column indices for the design
    col_start <- 2 + (d-1)*8
    col_end <- col_start + 7
    
    # Extract and process data
    bias_data <- result[[result_idx]][,col_start:col_end] - truepce[if(c == 1) 3 else 4]
    df_temp <- as.data.frame(bias_data)
    colnames(df_temp) <- c("PAR-0.5", "PAR-1", "PAR-2", "PAR-INF", 
                           "ML-0.5", "ML-1", "ML-2", "ML-INF")
    
    # Add metadata columns
    df_temp$design <- designs[d]
    df_temp$condition <- conditions[c]
    df_temp$row_id <- 1:nrow(df_temp)
    
    # Determine which estimators should be orange
    if (conditions[c] == "Always-takers") {
      if (designs[d] %in% c("(i)", "(ii)", "(iii)")) {
        highlight_est <- c("PAR-INF", "ML-INF")
      } else {
        highlight_est <- c("ML-INF")
      }
    } else {  # Compliers
      if (designs[d] %in% c("(i)", "(ii)", "(iii)")) {
        highlight_est <- c("PAR-INF", "ML-INF")
      } else {
        highlight_est <- c("ML-INF")
      }
    }
    
    # Reshape to long format
    df_long <- reshape2::melt(df_temp, 
                              id.vars = c("row_id", "design", "condition"),
                              variable.name = "estimator",
                              value.name = "bias")
    
    # Add highlight column
    df_long$highlight <- df_long$estimator %in% highlight_est
    
    combined_df <- rbind(combined_df, df_long)
  }
}

# save pdf
pdf("boxINF_00.pdf", width = 5, height = 8)  # Specify the size here
# Create the plot
ggplot(combined_df[which(combined_df$condition=="00"),], aes(x = estimator, y = bias, fill = highlight)) +
  geom_boxplot() +
  facet_grid(design ~ condition) +
  coord_cartesian(ylim = c(-1.6, 1.6)) +
  scale_fill_manual(values = c("gray", "orange")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "lightblue", color = "black"),  # Set background color
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    axis.title = element_text(size = 14)
  ) +
  labs(x = "", y = "Bias")
dev.off()  # Close the PDF device
