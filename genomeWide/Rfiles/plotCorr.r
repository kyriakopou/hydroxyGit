#read the .mat files
path <- "~/Desktop/"
fileName <- file.path(path, "Corr_d0_WT.csv")
Corr_d0 <- read.csv(fileName, header = FALSE)

#transform lists into matrices with names
vars <- c("uu", "um", "toth", "mm", "maint", 'deNovo', "hydroxy");

Corr_d0 <- matrix(unlist(Corr_d0), ncol = 7, byrow = TRUE, dimnames = list(vars, vars))
# Corr_d3 <- matrix(unlist(Corr_d3), ncol = 7, byrow = TRUE, dimnames = list(vars, vars))
# Corr_d6 <- matrix(unlist(Corr_d6), ncol = 7, byrow = TRUE, dimnames = list(vars, vars))

# Corr_d0_TKO <- matrix(unlist(Corr_d0), ncol = 7, byrow = TRUE, dimnames = list(vars, vars))
# Corr_d4_TKO <- matrix(unlist(Corr_d4), ncol = 7, byrow = TRUE, dimnames = list(vars, vars))
# Corr_d7_TKO <- matrix(unlist(Corr_d7), ncol = 7, byrow = TRUE, dimnames = list(vars, vars))
corrplot(Corr_d0, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
         )


