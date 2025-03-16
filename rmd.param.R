library(rmarkdown)

# Define the list of individuals
individuals <- data.frame(
  name = c("Soumyajit Roy", "Sudipta Sarkar","Rupanjan Mukherjee","Arnab Rakshit","Rhitankar Bandyopadhyay","Swarnadeep Datta","Subhrangsu Bhunia", "Arpan Dutta"),
  roll = c("MD2319", "MD2325","MD2314","MD2303","MD2313","MD2326","MD2324","MD2305"),
  email = c("soumyajitroy25356@gmail.com", "sudiptasarkar8607@gmail.com","rup629063@gmail.com","arnabrakshit2018@gmail.com","rhitankarbandyopadhyay@gmail.com","swarnadeepdatta0122@gmail.com","bhuniasubhrangsu19@gmail.com", "arpandutta0405@gmail.com")
)

# Path to  Rmd file
rmd_file <- "C:/Users/HP/OneDrive/Desktop/ComputingII/Assignment.param.Rmd"  

# Loop through each individual and generate a PDF
for (i in 1:nrow(individuals)) {
  output_file <- paste0("Density_Estimation_", individuals$roll[i], ".pdf")
  
  render(
    input = rmd_file,
    output_file = output_file,
    params = list(
      name = individuals$name[i],
      roll = individuals$roll[i],
      email = individuals$email[i]
    )
  )
  
  cat(sprintf("Generated PDF for %s: %s\n", individuals$name[i], output_file))
}






