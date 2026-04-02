
#' Function for diagnosing errors in a run of the pipeline
#'
#' @param input A character path to a project folder or a gpa object
#'
#' @importFrom dplyr bind_rows
#' @importFrom cli cli_abort cli_warn cli_inform cli_alert no qty
#' @export
#'
#' @author Zach Vig
diagnose_pipeline <- function(input) {
  
  old_warn <- getOption("warn")
  if (old_warn < 1) {
    on.exit(options(warn = old_warn), add = TRUE)
    options(warn = 1)
  }
  
  cli::cli_inform(
    c("Running pipeline diagnosis...", 
      "(Hit ESC at any time to cancel)")
    )
  
  if (checkmate::test_class(input, "character")) {
    if (isFALSE(checkmate::test_directory_exists(input))) {
      cli::cli_abort(
        "Provided directory does not exist: {.path {input}}"
      )
    }
    
    # if input is path to project directory
    proj_dir <- input
    proj_files <- list.files(proj_dir, include.dirs = TRUE)
    scheduler_scripts_dir <- file.path(input, "scheduler_scripts")
    sqlite_db <- list.files(proj_dir, pattern = "\\.sqlite$", full.names = TRUE)
    if (length(sqlite_db) > 1) {
      cli::cli_inform(
        c("More than one SQLite database found in folder.",
          "Which one would you like to use for diagnosis?")
      )
      db_choice <- prompt_input(
        type = "integer",
        required = TRUE,
        lower = 1,
        upper = length(sqlite_db),
        prompt = "Enter an integer"
      )
      sqlite_db <- sqlite_db[db_choice]
    }
    config_file <- file.path(proj_dir, "project_config.json")
    
    input_is_gpa <- FALSE
    
  } else if (checkmate::test_class(input, "glm_pipeline_arguments")) {
    
    # if input is gpa object
    proj_dir <- input$output_directory
    proj_files <- list.files(proj_dir, include.dirs = TRUE)
    scheduler_scripts_dir <- input$output_locations$scheduler_scripts
    sqlite_db <- input$output_locations$sqlite_db
    config_file <- input$output_locations$project_config_json
    
    input_is_gpa <- TRUE
    
  } else {
    cli::cli_abort(
      "Input must be a gpa object or a character path to a project/output directory."
    )
  }
  
  if (checkmate::test_directory_exists(scheduler_scripts_dir)) {
    
    scheduler_scripts <- list.files(scheduler_scripts_dir)
    
    # make sure scheduler script directory has files
    if (isFALSE(length(scheduler_scripts) > 0)) {
      cli::cli_abort(
        c("Scheduler scripts directory is empty.",
          "i" = "Maybe the analysis hasn't been run yet?")
      )
    }
    
  } else {
    cli::cli_abort(
      "Scheduler scripts directory not found: {.path {scheduler_scripts_dir}}"
    )
    # TODO: if SQLite doesn't exist, stop or return limited info?
  }
  
  if (checkmate::test_file_exists(sqlite_db)) {
    
    # make sure job tracking table exists in database
    if (isFALSE(sqlite_table_exists(sqlite_db, "job_tracking"))) {
      cli::cli_abort(
        c("Job tracking table does not exist in SQLite database.",
        "i" = "Maybe the analysis hasn't been run yet?")
      )
    }
    
  } else {
    cli::cli_abort(
      "SQLite database not found: {.path {sqlite_db}}"
    )
  }
  
  config_res <- summarize_project_config(config_file) # print submission history or return error code
  
  if (length(config_res) == 1 && config_res == "file_dne") {
    cli::cli_warn(
      c("Project configuration file does not exist.",
        "Batch folders will be read manually.")
    )
    
    ss_files <- list.dirs(scheduler_scripts_dir, full.names = TRUE)
    batch_dirs <- ss_files[grep("batch_", ss_files)]; n_batch <- length(batch_dirs)
    if (n_batch == 0) {
      cli::cli_abort(
        c("Batch folders not found.",
          "i" = "Could they have been moved or deleted from the scheduler scripts directory")
      )
    }
    sequence_ids <- sapply(basename(batch_dirs), function(x) gsub("batch_", "", x))
    # get creation times
    tmp <- sapply(
      batch_dirs, 
      function(x) file.info(x)$ctime, 
      simplify = FALSE,
      USE.NAMES = FALSE
      )
    creation_times <- do.call("c", tmp)
    ord <- sort(as.numeric(creation_times), index.return = TRUE, decreasing = TRUE)$ix
    
    # reorder so that display order matches selection indices
    batch_dirs <- batch_dirs[ord]
    sequence_ids <- sequence_ids[ord]
    creation_times <- creation_times[ord]
    
    tags <- c("{.emph (Latest)}", rep("", length(batch_dirs) - 1))
    
    cli::cli_inform(
      c("The pipeline for this analysis has been run {n_batch} time{?s}.",
        "These are the batch folders in the scheduler scripts directory:")
    )
    
    # list batch directory names in order of creation time (newest first)
    cli::cli_ol(
      cli::col_cyan(
        glue::glue("{basename(batch_dirs)} (Run on {format(as.POSIXct(creation_times), '%D')} 
                   at {format(as.POSIXct(creation_times), '%I:%M:%S %p')}) {tags}")
        )
    )
    
    default <- 1
    
  } else if (length(config_res) == 1 && config_res == "no_history") {
    
    cli::cli_abort(
      c("Project configuration file has no submission history.",
        "i" = "Maybe the analysis hasn't been run yet?")
    )
    
  } else {
    sequence_ids <- config_res
    default <- length(sequence_ids)
  }
  
  cli::cli_inform(
    "Which pipeline run would you like to diagnose?"
  )
  
  ss <- prompt_input(
    prompt = "Enter an integer",
    default = default,
    type = "integer",
    lower = 1,
    upper = length(sequence_ids),
    len = 1
  )
  
  this_sequence_id <- sequence_ids[ss]
  
  if (length(config_res) == 1 && config_res == "file_dne") {
    batch_directory <- batch_dirs[ss]
  } else {
    batch_directory <- get_sequence_info(config_file = config_file, sequence_id = this_sequence_id)$batch_directory
  }
  
  if(!checkmate::test_directory_exists(batch_directory)) {
    cli::cli_warn(
      c("Batch directory for this pipeline run does not exist.",
        "i" = "Could it have been moved or deleted from the scheduler scripts folder?")
    )
    batch_directory_exists <- FALSE
  } else {
    batch_directory_exists <- TRUE
  }
  
  # find sequence in SQLite database & make sure job was tracked
  tracking_df <- get_tracked_job_status(sequence_id = this_sequence_id, sqlite_db = sqlite_db)
  retrieval_time <- Sys.time()
  if (isTRUE(nrow(tracking_df) == 0)) {
    cli::cli_abort(
      "This pipeline run was not tracked in the SQLite database."
    )
  }
  
  # convert tracking data.frame to data.tree object
  this_sequence_tree <- tracking_df_to_tree(tracking_df)[[this_sequence_id]]
  n_chains <- this_sequence_tree$count
  
  # get pipeline steps
  chains <- lapply(paste0("chain", 1:n_chains), function(chain) this_sequence_tree[[chain]]$children)
  
  cli::cli_inform(
    "The run you selected had {n_chains} independent chain{?s} of steps:"
  )
  
  for (i_chain in 1:n_chains) {
    
    if (n_chains > 1) {
      cli::cli_h3(
        "Chain {i_chain}"
      )
    }
    
    # start chain summary list
    cli::cli_ol(id = "chain_summary")
    
    # list of steps for this chain
    steps <- chains[[i_chain]]
    
    for (i_step in seq_along(steps)) {
      
      this_step <- steps[[i_step]]
      
      # step information for printing
      this_step_name <- this_step$name
      this_step_title <- get_step_title(this_step_name)
      this_step_status <- this_step$status
      this_step_symbol <- get_status_symbol(this_step_status)
    
      if (this_step_status == "FAILED_BY_EXT") {
        # get upstream job name
        v <- sapply(steps, function(x) x$id == this_step$parent_id)
        upstream <- names(which(v))
        # status message to print
        status_message <- sprintf("FAILED (because `%s` failed)", get_step_title(upstream))
      } else {
        status_message <- this_step_status
      }
      
      n_children <- this_step$count
      bullet <- ifelse(isTRUE(n_children > 0), "\u00a6", "\u2218")
      
      # print step name and status
      cli::cli_li(
        "{this_step_title} [{this_step_symbol}] {status_message}"
      )
    
      # if step has children, print table of their statuses
      if (n_children > 0) {

        # get statuses
        children_status <- sapply(this_step$children, function(x) unname(x$Get("status")), simplify = TRUE)
        
        # format table of children statuses
        tab <- table(children_status)
        tab_codes <- names(tab)
        tab_labels <- sapply(tab_codes, function(x) ifelse(x == "FAILED_BY_EXT", "FAILED (due to parent)", x))
        n_row_tab <- length(tab)
        tab_count <- unname(tab)
        bullets <- c(rep("\u251c" , n_row_tab - 1), "\u2514")
        drop_rows <- sprintf("%s\u2500 %s of %s [%s] %s", bullets, tab_count, n_children, get_status_symbol(tab_codes), tab_labels)
        names(drop_rows) <- " "
        
        # print results
        cli::cli_bullets(
          c(" " = "Child Jobs:", drop_rows)
        )
        
      } 
      
    }
    
    cli::cli_end(id = "chain_summary")
    
  }
  
  cli::cli_inform(
    c("Would you like to examine any of these steps more closely?",
      "i" = "Type 'yes' to continue, or 'no' to return full job tree and exit now.")
  )
  continue <- prompt_input(
    type = "flag",
    required = TRUE
  )
  
  if (isFALSE(continue)) {
    # don't continue sequence if requested
    cli::cli_inform(
      c("Exiting diagnosis and returning full tracking tree...",
        "i" = "Make sure you load the {.pkg data.tree} package")
    )
    return(this_sequence_tree)
  }
  
  # otherwise, continue sequence
  if (n_chains > 1) {
    cli::cli_inform(
      c("Which chain would you like to examine more closely?",
      "i" = "There are {n_chains} chains")
    )
    chain_choice <- prompt_input(
      type = "integer",
      required = TRUE,
      lower = 1,
      upper = n_chains,
      prompt = "Enter an integer"
    )
  } else {
    chain_choice <- 1
  }
  
  n_steps <-length(chains[[chain_choice]])
  
  if (n_steps > 1) {
    cli::cli_inform(
      c("Which step would you like to examine more closely?",
      "i" = "There are {n_steps} steps in chain {chain_choice}")
    )
    step_choice <- prompt_input(
      type = "integer",
      required = TRUE,
      lower = 1,
      upper = n_steps,
      prompt = "Enter an integer"
    )
  } else {
    cli::cli_inform(
      c("There is only one step in this chain.",
        "Defaulting to that step...")
    )
    step_choice <- 1
  }
  
  this_step <- chains[[chain_choice]][[step_choice]]
  
  # this step info
  this_step_name <- this_step$name
  this_step_title <- get_step_title(this_step_name)
  this_step_status <- this_step$status
  
  if (this_step_status == "COMPLETED") {
    status_list <- c(
      "...was submitted at {this_step$time_submitted}",
      "...was started at {this_step$time_started}",
      "...{.emph successfully completed} at {this_step$time_ended}"
    )
  } else if (this_step_status == "QUEUED") {
    status_list <- c(
      "...was submitted at {this_step$time_submitted}",
      "...has {.emph not started} as of {retrieval_time}"
    )
  } else if (this_step_status == "STARTED") {
    status_list <- c(
      "...was submitted at {this_step$time_submitted}",
      "...was started at {this_step$time_started}",
      "...has {.emph not finished} as of {retrieval_time}"
    )
  } else {
    
    if (this_step_status == "FAILED_BY_EXT") {
      # get upstream job
      v <- sapply(steps, function(x) x$id == this_step$parent_id)
      upstream <- names(which(v))
      
      extra <- "because upstream job {.code {upstream}} failed"
    } else {
      extra <- ""
    }
    
    was_queued <- !is.na(this_step$time_submitted)
    was_started <- !is.na(this_step$time_started)
    
    if (isTRUE(was_queued)) {
      if (isTRUE(was_started)) {
        status_list <- c(
          "...was submitted at {this_step$time_submitted}",
          "...and was started at {this_step$time_started}",
          paste("...{.emph failed} before finishing", extra)
        )
      } else {
        status_list <- c(
          "...was submitted at {this_step$time_submitted}",
          paste("...{.emph failed} before starting", extra)
        )
      }
    } else {
      status_list <- c(
        paste("...{.emph failed} before being submitted", extra)
      )
    }

  }
  
  # downstream summary
  v <- sapply(steps, function(x) isTRUE(x$parent_id == this_step$id))
  downstream <- names(which(v)); n_downstream <- length(downstream) 
  
  downstream_message <- c(
    "...had {no(n_downstream)} downstream (dependent) job{?s/:/s:}"
  )
  # children summary
  children_message <- c(
    "...had {no(this_step$count)} child job{?s}"
  )
  
  # print this step detailed status message
  cli::cli_h3(
    "Job {.code {this_step_name}}..."
  )
  cli::cli_ul(id = "step_summary")
  cli::cli_li(
    c(status_list,
    downstream_message)
  )
  ds <- cli::cli_ul(); cli::cli_li(sprintf("{.code %s}", downstream)); cli::cli_end(ds)
  cli::cli_li(
    children_message,
    id = "step_summary"
  )
  cli::cli_end("step_summary")
  
  # initialize action choice set
  actions <- c("exit" = cli::col_cyan("Exit diagnosis and return full tracking tree"))
  max_choice <- 1L

  # get .Rout file for this step
  Rout <- paste0(this_step$compute_file, "out")
  
  # check if .Rout file exists
  if (checkmate::test_file_exists(Rout)) {
    Rout_exists <- TRUE
  } else  if (!checkmate::test_file_exists(Rout) && isTRUE(batch_directory_exists)) {
    tryRout <- file.path(batch_directory, paste0("batch_run_", this_step_name, ".Rout"))
    if (checkmate::test_file_exists(tryRout)) {
      Rout <- tryRout
      Rout_exists <- TRUE
    } else {
      Rout_exists <- FALSE
    }
  } else {
    Rout_exists <- FALSE
  }
  
  if (isFALSE(Rout_exists)) {
    cli::cli_warn(
      "Could not find this step's output file.",
      "i" = "You will not be able to view it."
    )
  } else {
    actions <- append(
      actions,
      c("view_Rout" = cli::col_cyan("View output file in console"),
        "return_Rout" = cli::col_cyan("Return output file as character object")),
      after = 0
    )
    max_choice = max_choice + 1
  }
  
  # add view children action?
  if (this_step$count > 0) {
    actions <- append(
      actions,
      c("view_children" = cli::col_cyan("View child job details (statuses, etc.)")),
      after = length(actions) - 1
    ) 
    max_choice <- max_choice + 1
  }
  
  # get feat table and log, if applicable
  if (isTRUE(this_step_name %in% c("setup_l1", "setup_run_l2", "setup_run_l3"))) {
    
    this_step_model_no <- readr::parse_number(this_step_name)
    
    if (isTRUE(input_is_gpa)) {
      fsl_table <- input[[glue::glue("l{this_step_model_no}_model_setup")]]$fsl
      fsl_table_exists <- TRUE
    } else {
      tryRData <- file.path(batch_directory, "run_pipeline_cache.RData")
      if (checkmate::test_file_exists(tryRData)) {
        tryCatch({
          gpa_obj <- local({load(tryRData); environment()})$gpa
          fsl_table <- gpa_obj[[glue::glue("l{this_step_model_no}_model_setup")]]$fsl
          fsl_table_exists <- TRUE
        },
        error = function (e) {
          cli::cli_warn(
            c("Unable to load the cached gpa object in batch directory.",
              "i" = "You will not be able to view feat status data.")
          )
          fsl_table_exists <- FALSE
        })
        
      } else {
        cli::cli_warn(
          c("Could not find cached gpa object in batch directory.",
            "i" = "You will not be able to view feat status data.")
        )
        fsl_table_exists <- FALSE
      }
    }
    
    if (isTRUE(input_is_gpa)) {
      try_fsl_log <- input$output_locations[[glue::glue("setup_l{this_step_model_no}_log_txt")]]
    } else {
      try_fsl_log <- file.path(proj_dir, glue::glue("logs/setup_l{this_step_model_no}_models.txt"))
    }
    
    if (checkmate::test_file_exists(try_fsl_log)) {
      fsl_log_exists <- TRUE
    } else {
      cli::cli_warn(
        c("Could not find level {this_step_model_no} models setup log.",
          "i" = "You will not be able to view it.")
      )
      fsl_log_exists <- FALSE
    }
    
  } else {
    fsl_log_exists <- FALSE
    fsl_table_exists <- FALSE
  }
  
  # update action set
  if (isTRUE(fsl_table_exists)) {
    actions <- append(
      actions,
      c("view_fsl_table" = cli::col_cyan("View feat status table")),
      after = length(actions) - 1
    ) 
    max_choice <- max_choice + 1
    
    total_feat <- nrow(fsl_table)
    n_feat_complete <- sum(fsl_table$feat_complete)
    n_feat_failed <- sum(fsl_table$feat_failed)
  }
  if (isTRUE(fsl_log_exists)) {
    actions <- append(
      actions,
      c("view_fsl_log" = cli::col_cyan("View feat log")),
      after = length(actions) - 1
    ) 
    max_choice <- max_choice + 1
  }
  
  # exit if no choices are available
  if (max_choice == 1L) {
    cli::cli_inform(
      c("Exiting diagnosis since no actions are available for this step...",
        "i" = "Re-run the function to diagnose another step")
    )
    return(invisible(NULL))
  }
  
  
  # start action loop (for now, never ends)
  action_loop <- TRUE
  
  while (isTRUE(action_loop)) {
    
    cli::cli_inform(
      c("Which action would you like to perform on job {.code {this_step_name}}?",
        "(Hit ESC to exit + Re-run function to examine a different step)"),
      wrap = TRUE
    )
    
    # display action choices
    cli::cli_ol(
      actions
    )
    
    action_choice <- prompt_input(
      prompt = "Enter an integer or hit ESC",
      type = "integer",
      lower = 1,
      upper = max_choice,
      default = 1,
      required = TRUE
    )
    
    if (names(actions)[action_choice] == "exit") {
      
      cli::cli_inform(
        c("Exiting diagnosis and returning full tracking tree...",
        "i" = "Make sure you load the {.pkg data.tree} package")
      )
      return(this_sequence_tree)
      
    } else if (names(actions)[action_choice] == "view_children") {
      # get children statuses
      children_status <- sapply(this_step$children, function(x) unname(x$Get("status")), simplify = TRUE)
      
      # format table of children statuses
      tab <- table(children_status)
      tab_count <- unname(tab)
      tab_name <- names(tab)
      tab_tag <- vapply(
        tab_name,
        function(x) switch(
          x,
          "COMPLETED"     = "{?has/have} {.emph completed successfully}",
          "STARTED"       = "{?has/have} {.emph started} by not finished:",
          "QUEUED"        = "{?has/have} been {.emph submitted} but have not started:",
          "FAILED"        = "{?has/have} {.emph failed}:",
          "FAILED_BY_EXT" = "{?has/have} {.emph failed} because the parent job failed:"
        ),
        character(1L)
      )

      # print child job statuses
      cli::cli_h3(
        "Job {.code {this_step_name}} child jobs:"
      )
      for (i in seq_along(tab_tag)) {
        
        j <- unname(which(children_status == tab_name[i]))
        jobs <- names(children_status[j])
        statuses <- unname(children_status[j])
        
        cli::cli_bullets(
          paste("{tab_count[i]} of {this_step$count} {qty(tab_count[i])}", tab_tag[i])
        )
        cli::cli_ol(id = "children_list")
        cli::cli_li(
          sprintf("{.code %s} [%s]", jobs, get_status_symbol(statuses))
        )
        cli::cli_end(id = "children_list")
      
      }
      
    } else if (names(actions)[action_choice] == "view_fsl_table") {
      
      cli::cli_h3(
        "Feat status for level {this_step_model_no}:"
      )
      cli::cli_ul(id = "feat_status")
      cli::cli_li(
        c("{n_feat_complete} of {total_feat} [{get_status_symbol('COMPLETED')}] COMPLETED ",
          "{n_feat_failed} of {total_feat} [{get_status_symbol('FAILED')}] FAILED ")
        
      )
      if (n_feat_failed > 0) {
        which_feat_failed <- which(fsl_table$feat_failed)
        feat_failed_ids <- fsl_table[which_feat_failed, c("id", "session", "run_number")]
        feat_failed_print <- apply(feat_failed_ids, 1, 
                                   function(row) glue::glue("ID {row[1]}, Session {row[2]}, Run {row[3]}"))
        cli::cli_ul(id = "failed_list")
        cli::cli_li(
          feat_failed_print
        )
        cli::cli_end(id = "failed_list")
      }
      cli::cli_end(id = "feat_status")
      
      
    } else if (names(actions)[action_choice] == "view_fsl_log") {
      
      cli::cli_inform(
        c("How many lines (from the end of the file) do you wish to view?")
      )
      
      lines_choice <- prompt_input(
        prompt = "Enter an integer",
        type = "integer",
        min = 1,
        default = 15,
        required = TRUE
      )
      
      view_log(input, from_end = TRUE, lines = lines_choice, level = this_step_model_no)
      
      continue <- TRUE
      
      while (isTRUE(continue)) {
        
        # view individual-level log
        cli::cli_inform(
          "Would you like to view an individual-level log?"
        )
        continue <- prompt_input(
          type = "flag",
          required = TRUE
        )
        if (isTRUE(continue)) {
          cli::cli_inform(
            "What is the ID of the log you would like to view?"
          )
          id <- prompt_input(
            prompt = "Enter an integer",
            type = "integer",
            len = 1
          )
          cli::cli_inform(
            c("How many lines (from the end of the file) do you wish to view?")
          )
          lines_choice <- prompt_input(
            prompt = "Enter an integer",
            type = "integer",
            min = 1,
            default = 15,
            required = TRUE
          )
          view_log(input, from_end = TRUE, lines = lines_choice, level = this_step_model_no, id = id)
        }
        
      }
      
    } else {
        
      # read in .Rout file and substitute warnings() instances (they cause problems)
      Rout_lines <- gsub("warnings\\(\\)", "`warnings`", readLines(Rout))
      Rout_n_lines <- length(Rout_lines)
      
      if (names(actions)[action_choice] == "return_Rout") {
        cli::cli_inform(
          "Exiting diagnosis and returning output file for job {.code {this_step_name}}..."
        )
        return(Rout_lines)
      } else {
        
        cli::cli_inform(
          c("How many lines (from the end of the file) do you wish to view?")
        )
        
        lines_choice <- prompt_input(
          prompt = "Enter an integer",
          type = "integer",
          lower = 1,
          default = min(Rout_n_lines, 15),
          required = TRUE
        )
        
        if (lines_choice > Rout_n_lines) {
          lines_choice <- Rout_n_lines
        }
        
        cli::cli_h1("Output file for job {.code {this_step_name}}")
        cli::cli_verbatim(
          c("...", tail(Rout_lines, lines_choice))
        )
        cli::cli_h1("")
        if (isFALSE(lines_choice == Rout_n_lines)) {
          cli::cli_bullets(
            c("i" = "{.emph Only showing the last {lines_choice} line{?s} of {Rout_n_lines} total lines}")
          )

        }
        
      }
      
    }
    
  }

}

#' Function for quick viewing of log files in the command line
#'
#' @param input A character path to a project folder or a gpa object.
#' @param from_end Logical. If TRUE, view lines from end of file instead of beginning.
#' @param lines An integer. How many lines from beginnning/end to view.
#' @param level An integer equal to 1, 2 or 3 (model level).
#' @param id An integer or NULL. The ID number if wishing to view individual-level log.
#'
#' @returns invisible NULL
#' @export
#'
#' @author Zach Vig
view_log <- function(input, from_end = TRUE, lines = 15L, level = NULL, id = NULL) {
  if (checkmate::test_class(input, "character")) {
    if (isFALSE(checkmate::test_directory_exists(input))) {
      cli::cli_abort(
        "Provided directory does not exist: {.path {input}}"
      )
    }
    
    # if input is path to project directory
    proj_dir <- input
    proj_files <- list.files(proj_dir, include.dirs = TRUE)
    input_is_gpa <- FALSE
    
  } else if (checkmate::test_class(input, "glm_pipeline_arguments")) {
    
    # if input is gpa object
    proj_dir <- input$output_directory
    input_is_gpa <- TRUE
    
  } else {
    cli::cli_abort(
      "Input must be a gpa object or a character path to a project/output directory."
    )
  }
  
  if(!checkmate::test_true(level %in% 1:3)){
    cli::cli_abort(
      "{.code level} must be 1, 2, or 3."
    )
  }
  
  if (!is.null(id)) {
    tag <- ", ID {id}"
  } else {
    tag <- ""
  }
  
  # get log directory
  if (isTRUE(input_is_gpa)) {
    log_dir <- input$output_locations$log_directory
  } else {
    log_dir <- file.path(input, "logs")
  }
  if (isFALSE(checkmate::test_directory_exists(log_dir))) {
    cli::cli_warn(
      c("Could not find log directory.",
        "i" = "You will not be able to view logs.")
    )
    return(invisible(FALSE))
  }
  
  # get log path
  if (!is.null(id)) {
    log_path <- file.path(log_dir, glue("subj{id}/setup_l{level}_models_subj{id}.txt"))
  } else {
    log_path <- file.path(log_dir, glue("setup_l{level}_models.txt"))
  }
  if (!file.exists(log_path)) {
    cli::cli_warn(
      c("Could not find level {level} models{glue(tag)} setup log.",
        "i" = "You will not be able to view it.")
    )
    return(invisible(FALSE))
  }

  # read in log file
  log_lines <- readLines(log_path)
  log_n_lines <- min(lines, length(log_lines)) # make sure lines choice is not too long
  
  # print lines
  cli::cli_h1("Log file for level {level} models{glue(tag)}")
  if (isTRUE(from_end)) {
    cli::cli_verbatim(
      c("...", tail(log_lines, log_n_lines))
    )
  } else {
    cli::cli_verbatim(
      c("...", head(log_lines, log_n_lines))
    )
  }
  cli::cli_h1("")
  if (isFALSE(lines == log_n_lines)) {
    cli::cli_bullets(
      c("i" = "{.emph Only showing the last {lines} line{?s} of {log_n_lines} total lines}")
    )
  }
  return(invisible(TRUE))
}

#' helper for printing a more formal title for the steps in the pipeline
#'
#' @param step Character string step name (as defined in run_glm_pipeline), e.g., 'run_l1'
#'
#' @keywords internal
get_step_title <- Vectorize(
  function(step) {
    switch(step, 
      "finalize_configuration" = "Finalize Configuration",
      "setup_l1" = "Set Up L1 Models",
      "split_backend_caches" = "Split Backend Caches",
      "run_l1_fsl" = "Run L1 Models",
      "setup_run_l2" = "Set Up/Run L2 Models",
      "setup_run_l3" = "Set Up/Run L3 Models",
      "cleanup_glm" = "Clean Up Pipeline",
      paste0("`", step, "`"))
  },
  USE.NAMES = FALSE
  )

#' helper for printing symbols based on status
#'
#' @param status Character string job status
#' @importFrom cli col_green col_magenta col_yellow col_red
#'
#' @keywords internal
get_status_symbol <- Vectorize(
  function(status) {
    switch(status,
           "COMPLETED" = cli::col_green("\u2714"), # checkmark
           "QUEUED" = cli::col_magenta("\u2197"), # arrow
           "STARTED" = cli::col_yellow("\u22ef"), # ellipsis
           "FAILED" = cli::col_red("\u2717"), # X mark
           "FAILED_BY_EXT" = cli::col_red("\u2717")) # X mark
  },
  USE.NAMES = FALSE
  )

