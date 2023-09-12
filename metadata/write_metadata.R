### quick script for writing SHAI metadata in R to a .txt file

data.pred.SHMI <- readRDS("metadata/data.pred.SHMI.RDS")

write.table(data.pred.SHMI, file = "metadata/data.pred.SHMI.txt", sep = "\t",
            row.names = FALSE)