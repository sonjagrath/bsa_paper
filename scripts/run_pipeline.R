##############################################
# run_pipeline.R
# Master pipeline to run the full QTL-seq analysis workflow
# Automatically logs output and saves key R objects
##############################################

# ---- 1. Load helpers ----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

# ---- 2. Define paths ----
scripts_dir <- "scripts"
results_dir <- "results"
data_dir <- "data"
objects_dir <- "objects"
log_dir <- file.path(results_dir, "logs")

dirs <- c(scripts_dir, results_dir, data_dir, objects_dir, log_dir)
for (dir in dirs) if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

# ---- 3. Start logging ----
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_dir, paste0("pipeline_log_", timestamp, ".txt"))
log_con <- file(log_file, open = "a")
sink(log_con, split = TRUE)
sink(log_con, type = "message", append = TRUE)

message("==================================================")
message("🚀 Starting QTL-seq analysis pipeline")
message(paste("Timestamp:", timestamp))
message("🖥 Session info:")
message(capture.output(sessionInfo()))
message("==================================================\n")

# ---- 4. Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
step <- NULL
if (length(args) > 0 && "--step" %in% args) {
  idx <- which(args == "--step") + 1
  if (idx <= length(args)) step <- as.integer(args[idx])
}

# ---- 5. Helper function to run scripts ----
run_script <- function(script_name, scripts_dir = "scripts") {
  message("\n▶ Running:", script_name)
  env <- new.env()
  script_path <- file.path(scripts_dir, script_name)
  
  if (!file.exists(script_path)) return(list(status = "failed", error = "Script not found"))
  
  result <- tryCatch({
    source(script_path, local = env)
    message("✅ Completed:", script_name)
    
    list(status = "success", error = NULL)
  },
  error = function(e) {
    message(paste("❌ Error in", script_name, ":", e$message))
    list(status = "failed", error = e$message)
  })
  
  return(result)
}

# ---- 6. Define workflow steps ----
workflow_steps <- list(
  "1" = list(script = "01_QTLseq_Analysis.R"),
  "2" = list(script = "02_QTL_Annotation_and_GO_Enrichment.R"),
  "3" = list(script = "03_Supplementary_Tables_and_Reports.R"))

# ---- 7. Run step(s) ----
start_time <- Sys.time()
results_summary <- list()

if (is.null(step)) {
  message("Running full pipeline (steps 1–3)...")
  for (s in names(workflow_steps)) {
    step_info <- workflow_steps[[s]]
    results_summary[[step_info$script]] <- run_script(step_info$script)
  }
} else if (as.character(step) %in% names(workflow_steps)) {
  step_info <- workflow_steps[[as.character(step)]]
  results_summary[[step_info$script]] <- run_script(step_info$script)
} else {
  message("⚠️ Invalid step argument. Valid steps: 1, 2, 3")
}

end_time <- Sys.time()
message("\n✅ Pipeline finished")
message(paste("Elapsed time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

sink(type = "message")
sink()
close(log_con)

###EOF
