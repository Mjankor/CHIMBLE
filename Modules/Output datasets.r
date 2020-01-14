
#************************** Final Output *********************************


#Output data to excel files.
Output <- function(){

if (exists("ScenarioPath")){
message ("Generating output files for ",  substr(ScenarioFile,1,nchar(ScenarioFile)-2))
write.table(Met, paste ("/Scenarios/", substr(ScenarioFile,1,nchar(ScenarioFile)-2), "/", Met_OutputPath, ".txt", sep = ""), sep="\t", row.names=FALSE, quote = FALSE)
write.table(Lake, paste ("/Scenarios/", substr(ScenarioFile,1,nchar(ScenarioFile)-2), "/", Lake_OutputPath, ".txt", sep = ""), sep="\t", row.names=FALSE, quote = FALSE)
write.table(Clean, paste ("/Scenarios/", substr(ScenarioFile,1,nchar(ScenarioFile)-2), "/", Clean_OutputPath, ".txt", sep = ""), sep="\t", row.names=FALSE, quote = FALSE)
write.table(Isotopes, paste ("/Scenarios/", substr(ScenarioFile,1,nchar(ScenarioFile)-2), "/", Isotopes_OutputPath, ".txt", sep = ""), sep="\t", row.names=FALSE, quote = FALSE)
write.table(Chemistry, paste ("/Scenarios/", substr(ScenarioFile,1,nchar(ScenarioFile)-2), "/", Chemistry_OutputPath, ".txt", sep = ""), sep="\t", row.names=FALSE, quote = FALSE)
write.table(Flux, paste ("/Scenarios/", substr(ScenarioFile,1,nchar(ScenarioFile)-2), "/", Flux_OutputPath, ".txt", sep = ""), sep="\t", row.names=FALSE, quote = FALSE)

} else {
write.table(Met, Met_OutputPath, sep="\t", row.names=FALSE, quote = FALSE)
write.table(Lake, Lake_OutputPath, sep="\t", row.names=FALSE, quote = FALSE)
write.table(Clean, Clean_OutputPath, sep="\t", row.names=FALSE, quote = FALSE)
write.table(Isotopes, Isotopes_OutputPath, sep="\t", row.names=FALSE, quote = FALSE)
write.table(Chemistry, Chemistry_OutputPath, sep="\t", row.names=FALSE, quote = FALSE)
write.table(Flux, Flux_OutputPath, sep="\t", row.names=FALSE, quote = FALSE)
}
}