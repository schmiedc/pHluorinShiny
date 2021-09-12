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
# DEPENDENCIES:
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2018-08-07
#
# ============================================================================
# saving of the result tables in a csv

writeToCsv <- function(outdir, resultname, table.signal, table.background, finalTable){
  
  write.csv(table.signal, file = file.path(outdir,paste0( resultname, "_RawSignal.csv" ), fsep = .Platform$file.sep))
  write.csv(table.background, file = file.path(outdir,paste0( resultname, "_RawBackground.csv" ), fsep = .Platform$file.sep))
  write.csv(finalTable, file = file.path(outdir,paste0( resultname, "_Mean.csv" ), fsep = .Platform$file.sep))
  
}
