# Phase-gated runner. This invocation intentionally stops after input audit.
code_dir <- dirname(normalizePath(sys.frame(1)$ofile))
source(file.path(code_dir, "01_environment_check.R"))
source(file.path(code_dir, "02_input_audit.R"))
message("PHASE 2 COMPLETE. Stop for mandatory human review before merge.")
