library(xlsx)
# ============================================================================
#
#  DESCRIPTION: Data analysis for Fíji pHlorin workflow
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut f r Molekulare Pharmakologie (FMP)
#               Cellular Imaging - Core facility
#               Campus Berlin-Buch
#               Robert-Roessle-Str. 10
#               13125 Berlin, Germany
#
#         BUGS:
#        NOTES: 
# DEPENDENCIES: xlsx: install.packages("xlsx")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2018-08-07
#
# ============================================================================
# saving of the result tables in a csv

writeToCsv <- function(outdir, resultname, table.signal, table.background, finalTable, tau){
  
  write.csv(table.signal, file = file.path(outdir,paste0( resultname, "_RawSignal.csv" ), fsep = .Platform$file.sep))
  write.csv(table.background, file = file.path(outdir,paste0( resultname, "_RawBackground.csv" ), fsep = .Platform$file.sep))
  write.csv(finalTable, file = file.path(outdir,paste0( resultname, "_Mean.csv" ), fsep = .Platform$file.sep))
  write.csv(tau, file = file.path(outdir,paste0( resultname, "_Tau.csv" ), fsep = .Platform$file.sep))
  
}

writeToXlsx <- function(outdir, resultname, table.signal, table.background, finalTable, tau){
  
  # saving of the result tables in a xlsx
  write.xlsx(finalTable, file = file.path(outdir,paste0("Mean_",resultname,".xlsx")), sheetName="Total")
  write.xlsx(tau, file = file.path(outdir,paste0("Tau_",resultname,".xlsx")), sheetName="tau")
  # save raw data
  # write.xlsx(table.signal, file = file.path(outdir,paste0("RawSignal_",resultname,".xlsx")), sheetName="Signal")
  # write.xlsx(table.background, file = file.path(outdir,paste0("RawBackground_",resultname,".xlsx")), sheetName="Background")
  
  # write each name into a different xlsx sheet
  namecount <- as.data.frame(table(finalTable$name))
  
  for (names in namecount$Var1){
    
    pername.sheet <- subset(finalTable, name == names)
    write.xlsx(pername.sheet, file = file.path(outdir,paste0("Mean_SplitPerSheet_",resultname,".xlsx")), append=TRUE, sheetName=names)
    
  }
  
}
