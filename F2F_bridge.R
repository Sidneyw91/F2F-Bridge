#step1: Prepare the environment
input = read.csv("/path to data")
#read in table, with test sequences as rows and features as columns
Title = input[,1]
input[,1] = NULL
install.packages("fmsb")
install.packages("glmnet")
path = "/path to output directory"

input$Wet.Lab.Score.1 = NULL
input$Wet.Lab.Score.1.1 = NULL
####################################################################################
#step2: Define F2F_bridge function
F2F_bridge = function(input,row, plot = TRUE){
  library(fmsb)
  colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) )
  colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9))
  score = 0

    input = as.data.frame(t(input))
  
  alpha = as.numeric(paste(input$V3))
  beta = as.numeric(paste(input[,row]))
  
  ncolumns <- ncol(input)
  ar <- function(alpha, beta, input) {
    for (i in 1:nrow(input))
      score = score + abs(alpha[[i]]-beta[[i]])
      score = score/ncol(input)
    
    return(score)
    
  }
  
  score = ar(alpha, beta, input)
  
  input = as.data.frame(t(input))
  score = round(score, digits = 2)
  
  
  newList <- list("% outside allowable region" = ar, "protein score" = score)
  
  
  if (plot)
    graph = radarchart(input[c(1:3,row),], axistype=6,pty = 32, seg=10, plty=1, vlcex =  0.5,
                       pcol=colors_border , pfcol=colors_in , plwd=2,
                       #custom the grid
                       cglcol="grey", cglty=1, axislabcol="red", cglwd=0.8,
                       #custom labels
                       title = Title[row],
                       centerzero = TRUE
                       
    )
  legend(x=.7, y=1.2, legend = score, bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.7, pt.cex=3)
  return(newList)
  return(graph)
}
##########################################################################################################

#step 3: Use F2F-Bridge on inputted and prepared data frame
protein_stats = c()
for (protein in (4:nrow(input))){
  protein_stats = c(protein_stats,(F2F_bridge(input, protein)))
  print("plots are generating")
  plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
  plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
  file.copy(from=plots.png.paths, to= path)
   }

#####################################################################################################
  
#step 4: Extract the scores of the test sequences
Optimal_Protein = function(output){
df <- data.frame(matrix(unlist(output), nrow = length(4:nrow(input)), byrow=T))
rownames(df) = Title[4:nrow(input)]
df$X1 = NULL
df$X2 = as.numeric(df$X2)
print(rownames(df)[which(df == min(df))])
df$labels = rownames(df)
df <- df[order(df$X2),c(1,2)]
df$labels = NULL
fileName = paste(path, 'sequence_scores.txt')
write.table(df, fileName)
return(df)
}
print(Optimal_Protein(protein_stats))
######################################################################################################

#Step 5: No database method
lasso_select = function(data, results){
  input_lr = data[-c(1:3),]
  train = results
  library(glmnet)
  c = glmnet(as.matrix(input_lr), train, standardize = FALSE, alpha = 1)
  cvfit = cv.glmnet(as.matrix(input_lr), train)
  plot(c)
  plot(cvfit)
  cvfit$lambda.min
  result =coef(cvfit, s = "lambda.min")
  return(result)
}

