# _save_objects.R
# Helper function to save R objects with timestamped filenames

save_object <- function(obj, name, dir = "objects") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_path <- file.path(dir, paste0(name, "_", timestamp, ".RDS"))
  
  saveRDS(obj, file = file_path)
  message(paste("💾 Saved object:", name, "to", file_path))
  
  return(file_path)
}

