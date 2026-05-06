#===========================================================
# R/qtl + ASMap Analysis for Fire Blight Resistance
# Cleaned workflow:
#   - QC / diagnostics
#   - automatic allele-switch detection/correction
#   - physical map (Mb) for visualization
#   - ASMap genetic map (when possible)
#   - marker regression, manual Kruskal-Wallis, interval mapping, MQM
#   - permutation-based QTL scanning
#   - phenotype class summaries for resistant / intermediate / susceptible
#   - primary QTL peak genotype-effect summaries and plots
#   - no readable map creation
#   - no noisy warnings from fallback logic or MQM failures
#===========================================================

suppressPackageStartupMessages({
  library(qtl)
  library(ASMap)
  library(tidyverse)
  library(cli)
})

set.seed(12345)

# -----------------------------
# File paths
# -----------------------------
geno_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/R5xRG_Genome_Studio_Carolina_08_2016.txt"
pheno_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/Phenotype_Genotype_Mapping.txt"
snp_callrate_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/SNP_CallRates.txt"
output_dir <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Outputs"

all_markers_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/All_Markers.tsv"
marker_9k_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/9K_SNP_positions_DH_genome.txt"
mr5_map_file <- "C:/Users/StandishJ/Downloads/FinalProject/Rqtl_Inputs/Apple-IM-F1-Mr5_cM-Map.tsv"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Parameters
# -----------------------------
callrate_cutoff <- 0.95
error_prob <- 0.005
perm_n <- 10000

qtl_window_cm <- 5
candidate_keep_top_n_if_none <- 20

use_asmap <- TRUE
require_asmap_success <- FALSE
asmap_BC_gen <- 0L
asmap_F_gen <- 2L
asmap_p_value <- 2
asmap_noMap_dist <- 20
asmap_noMap_size <- 0
asmap_miss_thresh <- 0.2
asmap_trace <- FALSE
asmap_merge_gap <- 5

# -----------------------------
# Helpers
# -----------------------------
run_step <- function(label, expr) {
  cli_h2(label)
  start <- Sys.time()
  result <- tryCatch(
    withCallingHandlers(
      eval.parent(substitute(expr)),
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) invokeRestart("muffleMessage")
    ),
    error = function(e) stop(e)
  )
  elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
  cli_alert_success(paste0(label, " finished in ", elapsed, " sec"))
  result
}

quiet_eval <- function(expr) {
  val <- NULL
  ok <- tryCatch({
    capture.output({
      suppressWarnings(suppressMessages(
        val <<- eval.parent(substitute(expr))
      ))
    })
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  val
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.\\.[0-9]+\\.\\.", ".", x)
  x <- gsub("\\s*\\[.*?\\]", "", x)
  x <- gsub("\\.", "", x)
  trimws(x)
}

keep_nonempty_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

chr_main_from_label <- function(x) {
  x <- as.character(x)
  suppressWarnings(as.integer(stringr::str_extract(x, "^\\d+")))
}

sort_chr_labels <- function(x) {
  x <- unique(as.character(x))
  x[order(chr_main_from_label(x), x, na.last = TRUE)]
}

safe_pdf <- function(file, width = 11, height = 8.5, expr) {
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

safe_write_capture <- function(expr, file) {
  out <- tryCatch(
    capture.output(suppressWarnings(expr)),
    error = function(e) c(paste0("ERROR: ", conditionMessage(e)))
  )
  writeLines(out, con = file)
}

is_f2_cross <- function(cross) {
  ct <- tryCatch(qtl::crosstype(cross), error = function(e) NA_character_)
  !is.na(ct) && tolower(ct) == "f2"
}

ensure_map_names <- function(cross) {
  if (is.null(cross) || !inherits(cross, "cross")) return(cross)
  for (chr in names(cross$geno)) {
    m <- cross$geno[[chr]]$map
    nm <- names(m)
    if (is.null(nm) || any(is.na(nm) | nm == "")) {
      data_nm <- colnames(cross$geno[[chr]]$data)
      if (length(data_nm) == length(m)) {
        names(m) <- data_nm
        cross$geno[[chr]]$map <- m
      }
    }
  }
  cross
}

pick_chr_name <- function(cross_obj, chr_main) {
  chr_main <- as.character(chr_main)
  chr_names <- names(cross_obj$geno)
  exact <- chr_names[chr_names == chr_main]
  if (length(exact) > 0) return(exact[1])
  pat <- paste0("^", chr_main, "(\\.|$)")
  m <- chr_names[grepl(pat, chr_names)]
  if (length(m) > 0) return(m[1])
  NA_character_
}

get_geno_by_marker <- function(cross_obj, marker_id) {
  for (chr in names(cross_obj$geno)) {
    mat <- tryCatch(pull.geno(cross_obj, chr = chr), error = function(e) NULL)
    if (!is.null(mat) && marker_id %in% colnames(mat)) {
      return(mat[, marker_id])
    }
  }
  NULL
}

build_map_table <- function(cross, marker_lookup, pos_label = "cM") {
  map_list <- pull.map(cross)
  
  out <- bind_rows(lapply(names(map_list), function(chr_label) {
    map_vec <- map_list[[chr_label]]
    marker_names <- names(map_vec)
    if (is.null(marker_names)) marker_names <- rep(NA_character_, length(map_vec))
    tibble(
      chr_label = as.character(chr_label),
      marker_id = as.character(marker_names),
      map_pos = as.numeric(map_vec)
    )
  }))
  
  out %>%
    left_join(marker_lookup %>% select(marker_id, chr_main), by = "marker_id") %>%
    mutate(chr_main = coalesce(chr_main, chr_main_from_label(chr_label))) %>%
    rename(!!pos_label := map_pos)
}

parse_marker_annotations <- function(all_markers_file, marker_9k_file) {
  all_markers <- read.delim(
    all_markers_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  )
  
  marker_annot_clean <- all_markers %>%
    transmute(
      unique_name = trimws(`Unique Name`),
      marker_name = trimws(`Marker Name`),
      organism = trimws(`Organism`),
      mapped_organism = trimws(`Mapped Organism`),
      marker_type = trimws(`Marker Type`),
      linkage_group = trimws(`Linkage Group`),
      map_start = suppressWarnings(as.numeric(`Map Position (Start)`)),
      map_stop = suppressWarnings(as.numeric(`Map Position (Stop)`)),
      genome_start = suppressWarnings(as.numeric(`Genome Position (Start)`)),
      genome_stop = suppressWarnings(as.numeric(`Genome Position (Stop)`)),
      trait = Trait,
      publications = Publications
    ) %>%
    mutate(join_name = if_else(!is.na(marker_name) & marker_name != "", marker_name, unique_name)) %>%
    filter(!is.na(join_name), join_name != "") %>%
    distinct(join_name, .keep_all = TRUE)
  
  marker_9k <- read.delim(
    marker_9k_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  ) %>%
    rename(
      marker_name = SNPName,
      chr_9k = Chromosome,
      pos_9k = Position
    ) %>%
    mutate(
      marker_name = trimws(marker_name),
      chr_num_9k = suppressWarnings(as.integer(gsub("^Chr", "", chr_9k))),
      pos_bp_9k = suppressWarnings(as.numeric(pos_9k)),
      pos_mb_9k = pos_bp_9k / 1e6
    ) %>%
    filter(
      !is.na(marker_name),
      !is.na(chr_num_9k),
      !is.na(pos_bp_9k),
      chr_9k != "multihit",
      pos_9k != "-"
    ) %>%
    distinct(marker_name, .keep_all = TRUE)
  
  marker_annot_clean %>%
    left_join(
      marker_9k %>% select(marker_name, chr_num_9k, pos_bp_9k, pos_mb_9k),
      by = "marker_name"
    ) %>%
    mutate(
      chr_final = coalesce(chr_num_9k, suppressWarnings(as.integer(stringr::str_extract(linkage_group, "\\d+")))),
      pos_final_bp = coalesce(pos_bp_9k, map_start),
      pos_final_mb = pos_final_bp / 1e6
    )
}

extract_qtl_peak_effects <- function(cross_obj, peaks_tbl, map_tbl, pheno_col) {
  if (is.null(peaks_tbl) || nrow(peaks_tbl) == 0) {
    return(list(effect = tibble(), quartile = tibble(), plot_df = tibble()))
  }
  
  pheno_vec <- cross_obj$pheno[[pheno_col]]
  effect_list <- list()
  quartile_list <- list()
  plot_df_list <- list()
  
  for (i in seq_len(nrow(peaks_tbl))) {
    peak_chr_main <- peaks_tbl$chr_main[i]
    peak_pos <- as.numeric(peaks_tbl$pos[i])
    peak_label <- peaks_tbl$peak_label[i]
    
    chr_map <- map_tbl %>% filter(chr_main == peak_chr_main)
    if (nrow(chr_map) == 0) next
    
    nearest_marker <- chr_map %>%
      mutate(dist_to_peak = abs(cM_final - peak_pos)) %>%
      arrange(dist_to_peak, marker_id) %>%
      slice(1)
    if (nrow(nearest_marker) == 0) next
    marker_id <- nearest_marker$marker_id[1]
    
    geno_vec <- get_geno_by_marker(cross_obj, marker_id)
    if (is.null(geno_vec)) next
    
    df <- tibble(
      sample_id = names(geno_vec),
      genotype_code = as.character(geno_vec),
      phenotype = pheno_vec
    ) %>%
      filter(!is.na(genotype_code), !is.na(phenotype))
    
    if (nrow(df) == 0) next
    
    eff <- df %>%
      group_by(genotype_code) %>%
      summarise(
        n = n(),
        mean_pheno = mean(phenotype),
        sd_pheno = sd(phenotype),
        median_pheno = median(phenotype),
        .groups = "drop"
      ) %>%
      mutate(
        peak_label = peak_label,
        peak_chr_main = peak_chr_main,
        peak_pos = peak_pos,
        marker_id = marker_id
      )
    
    q25 <- as.numeric(quantile(df$phenotype, 0.25, na.rm = TRUE))
    q75 <- as.numeric(quantile(df$phenotype, 0.75, na.rm = TRUE))
    
    quart <- df %>%
      mutate(
        phenotype_group = case_when(
          phenotype <= q25 ~ "Resistant_Q1",
          phenotype >= q75 ~ "Susceptible_Q4",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(phenotype_group)) %>%
      count(genotype_code, phenotype_group, name = "n") %>%
      tidyr::pivot_wider(names_from = phenotype_group, values_from = n, values_fill = 0) %>%
      mutate(
        peak_label = peak_label,
        peak_chr_main = peak_chr_main,
        peak_pos = peak_pos,
        marker_id = marker_id
      )
    
    effect_list[[length(effect_list) + 1]] <- eff
    quartile_list[[length(quartile_list) + 1]] <- quart
    plot_df_list[[length(plot_df_list) + 1]] <- df %>%
      mutate(
        peak_label = peak_label,
        peak_chr_main = peak_chr_main,
        peak_pos = peak_pos,
        marker_id = marker_id
      )
  }
  
  list(
    effect = bind_rows(effect_list),
    quartile = bind_rows(quartile_list),
    plot_df = bind_rows(plot_df_list)
  )
}

merge_asmap_to_base_chromosomes <- function(cross_obj) {
  chr_groups <- names(cross_obj$geno)
  chr_base <- chr_main_from_label(chr_groups)
  if (all(is.na(chr_base))) return(cross_obj)
  
  merge_list <- split(chr_groups, chr_base)
  merge_list <- merge_list[lengths(merge_list) > 1]
  if (length(merge_list) == 0) return(cross_obj)
  
  ASMap::mergeCross(cross_obj, merge = merge_list, gap = asmap_merge_gap)
}

build_marker_map_long <- function(cross_obj) {
  map_list <- pull.map(cross_obj)
  
  bind_rows(lapply(names(map_list), function(chr_label) {
    map_vec <- map_list[[chr_label]]
    tibble(
      chr = as.character(chr_label),
      marker_id = as.character(names(map_vec)),
      pos = as.numeric(map_vec)
    )
  })) %>%
    filter(!is.na(marker_id), nzchar(marker_id), !is.na(pos))
}

manual_kruskal_scan <- function(cross_obj, pheno_col) {
  geno_mat <- pull.geno(cross_obj)
  pheno_vec <- cross_obj$pheno[[pheno_col]]
  
  map_tbl <- build_marker_map_long(cross_obj)
  idx <- match(colnames(geno_mat), map_tbl$marker_id)
  map_tbl <- map_tbl[idx, , drop = FALSE]
  
  out <- vector("list", ncol(geno_mat))
  
  for (j in seq_len(ncol(geno_mat))) {
    g <- geno_mat[, j]
    y <- pheno_vec
    ok <- complete.cases(g, y)
    
    g <- as.factor(g[ok])
    y <- y[ok]
    
    if (length(y) < 2 || nlevels(g) < 2) next
    
    kt <- tryCatch(kruskal.test(y ~ g), error = function(e) NULL)
    if (is.null(kt)) next
    
    p <- as.numeric(kt$p.value)
    p <- max(p, .Machine$double.xmin)
    
    out[[j]] <- tibble(
      marker_id = colnames(geno_mat)[j],
      chr = as.character(map_tbl$chr[j]),
      pos = as.numeric(map_tbl$pos[j]),
      n = length(y),
      groups = nlevels(g),
      statistic = unname(as.numeric(kt$statistic)),
      p_value = p,
      kw_score = -log10(p)
    )
  }
  
  bind_rows(out) %>%
    arrange(chr_main_from_label(chr), chr, pos, marker_id)
}

align_asmap_map_to_f2 <- function(f2_cross, asmap_cross) {
  if (is.null(f2_cross) || is.null(asmap_cross)) return(NULL)
  if (!inherits(f2_cross, "cross") || !inherits(asmap_cross, "cross")) return(NULL)
  
  asmap_map <- pull.map(asmap_cross)
  out <- list()
  
  for (chr_label in names(asmap_map)) {
    chr_main <- chr_main_from_label(chr_label)
    f2_chr <- pick_chr_name(f2_cross, chr_main)
    if (is.na(f2_chr)) return(NULL)
    
    f2_markers <- colnames(f2_cross$geno[[f2_chr]]$data)
    map_vec <- asmap_map[[chr_label]]
    map_names <- names(map_vec)
    
    if (is.null(map_names) || all(is.na(map_names)) || any(!nzchar(map_names))) {
      if (length(map_vec) == length(f2_markers)) {
        names(map_vec) <- f2_markers
      } else {
        return(NULL)
      }
    } else {
      if (all(f2_markers %in% map_names)) {
        map_vec <- map_vec[f2_markers]
      } else if (length(map_vec) == length(f2_markers)) {
        names(map_vec) <- f2_markers
      } else {
        return(NULL)
      }
    }
    
    out[[f2_chr]] <- map_vec
  }
  
  if (length(out) == 0) return(NULL)
  out
}

apply_map_list_to_cross <- function(cross_obj, map_list) {
  if (is.null(cross_obj) || is.null(map_list)) return(cross_obj)
  if (!inherits(cross_obj, "cross")) return(cross_obj)
  
  for (chr in names(map_list)) {
    if (chr %in% names(cross_obj$geno)) {
      map_vec <- map_list[[chr]]
      if (length(map_vec) == ncol(cross_obj$geno[[chr]]$data)) {
        names(map_vec) <- colnames(cross_obj$geno[[chr]]$data)
        cross_obj$geno[[chr]]$map <- map_vec
      }
    }
  }
  ensure_map_names(cross_obj)
}

run_mqm_analysis <- function(cross_obj, pheno_col) {
  if (!is_f2_cross(cross_obj)) return(NULL)
  
  cofactor_try <- c(20, 15, 10, 8, 5)
  for (ncof in cofactor_try) {
    mqm_cofactors <- tryCatch(
      quiet_eval(
        mqmautocofactors(
          cross_obj,
          num = ncof,
          distance = 5,
          dominance = FALSE,
          plot = FALSE,
          verbose = FALSE
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(mqm_cofactors)) next
    if (length(mqm_cofactors) == 0) next
    
    mqm_result <- tryCatch(
      quiet_eval(
        mqmscan(
          cross_obj,
          cofactors = mqm_cofactors,
          pheno.col = pheno_col,
          model = "additive",
          forceML = FALSE,
          cofactor.significance = 0.02,
          em.iter = 1000,
          window.size = 25.0,
          step.size = 5.0,
          logtransform = FALSE,
          estimate.map = FALSE,
          plot = FALSE,
          verbose = FALSE,
          outputmarkers = TRUE,
          multicore = TRUE,
          batchsize = 10,
          n.clusters = 1,
          test.normality = FALSE,
          off.end = 0
        )
      ),
      error = function(e) NULL
    )
    
    if (!is.null(mqm_result)) {
      return(list(result = mqm_result, cofactors = mqm_cofactors, ncof = ncof))
    }
  }
  
  NULL
}

# -----------------------------
# Read genotype and phenotype files
# -----------------------------
geno_raw <- run_step("Reading genotype file", {
  suppressWarnings(read.delim(geno_file, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE))
})

pheno <- run_step("Reading phenotype file", {
  read.delim(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
})

pheno$CleanID <- clean_id(pheno$CleanID)

# -----------------------------
# Sample columns
# -----------------------------
aux_col <- match("Aux", names(geno_raw))
if (is.na(aux_col)) stop("Could not find 'Aux' column in genotype file.")

sample_cols_raw <- names(geno_raw)[5:(aux_col - 1)]
sample_cols_raw <- keep_nonempty_chr(sample_cols_raw)

sample_cols_clean <- clean_id(sample_cols_raw)
sample_map <- setNames(sample_cols_raw, sample_cols_clean)

# -----------------------------
# SNP QC filter
# -----------------------------
snp_callrate_df <- run_step("Reading SNP call rate file", {
  read.table(snp_callrate_file, header = TRUE, stringsAsFactors = FALSE)
})

snps_keep <- snp_callrate_df$Name[snp_callrate_df$CallRate >= callrate_cutoff]
snps_keep <- keep_nonempty_chr(snps_keep)

# -----------------------------
# Match samples
# -----------------------------
common_ids <- intersect(names(sample_map), pheno$CleanID)
common_ids <- keep_nonempty_chr(common_ids)

if (length(common_ids) == 0) {
  stop("No matching sample IDs found between genotype and phenotype files after cleaning.")
}

sample_cols_keep_raw <- unname(sample_map[common_ids])
sample_cols_keep_raw <- keep_nonempty_chr(sample_cols_keep_raw)

if (length(sample_cols_keep_raw) == 0) {
  stop("No genotype sample columns matched the phenotype file.")
}

# -----------------------------
# Filter genotype rows/columns
# -----------------------------
geno_filt <- geno_raw[
  geno_raw$Name %in% snps_keep,
  c("Index", "Name", "Chr", "Position", sample_cols_keep_raw),
  drop = FALSE
]

if (ncol(geno_filt) <= 4) stop("No sample columns left in genotype data after filtering.")

geno_filt$Chr <- suppressWarnings(as.integer(geno_filt$Chr))
geno_filt$Position <- suppressWarnings(as.numeric(geno_filt$Position))

valid_chr <- !is.na(geno_filt$Chr) & geno_filt$Chr %in% 1:17
valid_pos <- !is.na(geno_filt$Position) & geno_filt$Position > 0

dropped_markers <- geno_filt[!(valid_chr & valid_pos), c("Index", "Name", "Chr", "Position"), drop = FALSE]
write.csv(
  dropped_markers,
  file = file.path(output_dir, "Dropped_Markers_Invalid_Chr_or_Position.csv"),
  row.names = FALSE,
  quote = FALSE
)

geno_filt <- geno_filt[valid_chr & valid_pos, , drop = FALSE]

names(geno_filt)[5:ncol(geno_filt)] <- common_ids

geno_filt <- geno_filt %>%
  arrange(Chr, Position, Index) %>%
  mutate(marker_id = sprintf("SNP_%05d", as.integer(Index)))

if (any(is.na(geno_filt$marker_id) | !nzchar(geno_filt$marker_id))) {
  stop("Some marker_id values are missing after creation.")
}
if (anyDuplicated(geno_filt$marker_id)) {
  stop("Duplicate marker_id values detected. Check the SNP Index column.")
}

marker_lookup <- geno_filt %>%
  transmute(
    Index,
    marker_id,
    SNPName = Name,
    chr_main = Chr,
    Position
  )

write.csv(
  marker_lookup,
  file = file.path(output_dir, "Marker_ID_Lookup.csv"),
  row.names = FALSE,
  quote = FALSE
)

physical_map_tbl <- marker_lookup %>%
  transmute(
    chr_main = chr_main,
    marker_id = marker_id,
    phys_bp = Position,
    phys_mb = Position / 1e6
  )

write.csv(
  physical_map_tbl,
  file = file.path(output_dir, "Physical_Map_Table.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Filter phenotype table
# -----------------------------
pheno_filt <- pheno[
  match(common_ids, pheno$CleanID),
  c("CleanID", "MeanDiseaseIncidence"),
  drop = FALSE
]

pheno_filt$MeanDiseaseIncidence <- as.numeric(pheno_filt$MeanDiseaseIncidence)
pheno_filt <- pheno_filt[!is.na(pheno_filt$MeanDiseaseIncidence), , drop = FALSE]

common_ids <- pheno_filt$CleanID
common_ids <- keep_nonempty_chr(common_ids)

geno_filt <- geno_filt[, c("Index", "Name", "marker_id", "Chr", "Position", common_ids), drop = FALSE]

# -----------------------------
# Phenotype classing for reporting / plot clarity
# -----------------------------
q25 <- as.numeric(quantile(pheno_filt$MeanDiseaseIncidence, 0.25, na.rm = TRUE))
q50 <- as.numeric(quantile(pheno_filt$MeanDiseaseIncidence, 0.50, na.rm = TRUE))
q75 <- as.numeric(quantile(pheno_filt$MeanDiseaseIncidence, 0.75, na.rm = TRUE))

pheno_class_df <- pheno_filt %>%
  mutate(
    Phenotype_Class = case_when(
      MeanDiseaseIncidence <= q25 ~ "Resistant",
      MeanDiseaseIncidence >= q75 ~ "Susceptible",
      TRUE ~ "Intermediate"
    ),
    Phenotype_Class = factor(Phenotype_Class, levels = c("Resistant", "Intermediate", "Susceptible"))
  )

write.csv(
  pheno_class_df,
  file = file.path(output_dir, "Phenotype_Class_Assignment.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  pheno_class_df %>%
    count(Phenotype_Class, name = "n") %>%
    mutate(q25 = q25, q50 = q50, q75 = q75),
  file = file.path(output_dir, "Phenotype_Class_Summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Phenotype ranking plot
# -----------------------------
pheno_plot_df <- pheno %>%
  filter(CleanID %in% common_ids | CleanID %in% c("R5GType", "RoyalGalaGType")) %>%
  mutate(
    CleanID = as.character(CleanID),
    MeanDiseaseIncidence = as.numeric(MeanDiseaseIncidence),
    group = if_else(CleanID %in% c("R5GType", "RoyalGalaGType"), "Parent", "Progeny")
  ) %>%
  filter(!is.na(MeanDiseaseIncidence)) %>%
  arrange(MeanDiseaseIncidence) %>%
  mutate(label = factor(CleanID, levels = CleanID))

p_pheno <- ggplot(pheno_plot_df, aes(x = MeanDiseaseIncidence, y = reorder(label, MeanDiseaseIncidence))) +
  geom_col(aes(fill = group), width = 0.8, alpha = 0.9) +
  scale_fill_manual(values = c("Parent" = "firebrick", "Progeny" = "steelblue")) +
  labs(
    title = "Progeny ranked from most resistant to least resistant",
    x = "Mean disease incidence (0 = most resistant)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(output_dir, "Phenotype_Ranked_Barplot.png"),
  plot = p_pheno,
  width = 12,
  height = max(6, 0.18 * nrow(pheno_plot_df)),
  dpi = 200
)

p_pheno_class <- ggplot(pheno_class_df, aes(x = Phenotype_Class, y = MeanDiseaseIncidence)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  labs(
    title = "Phenotype classes based on quartiles",
    x = NULL,
    y = "Mean disease incidence"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "Phenotype_Class_Boxplot.png"),
  plot = p_pheno_class,
  width = 9,
  height = 6,
  dpi = 200
)

# -----------------------------
# Write R/qtl input files
# -----------------------------
rqtl_geno_file <- file.path(output_dir, "data.csv")
rqtl_pheno_file <- file.path(output_dir, "pheno.csv")

phe_df <- data.frame(
  ID = common_ids,
  MeanDiseaseIncidence = pheno_filt$MeanDiseaseIncidence,
  stringsAsFactors = FALSE
)
write.csv(phe_df, file = rqtl_pheno_file, row.names = FALSE, quote = FALSE)

geno_mat <- as.matrix(geno_filt[, common_ids, drop = FALSE])
geno_mat <- toupper(trimws(geno_mat))
geno_mat[!(geno_mat %in% c("AA", "AB", "BB"))] <- NA

geno_indiv <- as.data.frame(t(geno_mat), stringsAsFactors = FALSE)
geno_indiv <- geno_indiv[common_ids, , drop = FALSE]
colnames(geno_indiv) <- geno_filt$marker_id

position_mb <- geno_filt$Position / 1e6

geno_lines <- c(
  paste(c("ID", geno_filt$marker_id), collapse = ","),
  paste(c("", as.character(geno_filt$Chr)), collapse = ","),
  paste(c("", as.character(position_mb)), collapse = ","),
  vapply(common_ids, function(id) {
    paste(c(id, as.character(geno_indiv[id, ])), collapse = ",")
  }, character(1))
)

writeLines(geno_lines, con = rqtl_geno_file)

rqtl_geno_file <- normalizePath(rqtl_geno_file, mustWork = TRUE)
rqtl_pheno_file <- normalizePath(rqtl_pheno_file, mustWork = TRUE)

message("Writing genotype file to: ", rqtl_geno_file)
message("Writing phenotype file to: ", rqtl_pheno_file)

# -----------------------------
# Load into R/qtl as the physical-position cross
# -----------------------------
cross_phys <- suppressWarnings(read.cross(
  format = "csv",
  file = rqtl_geno_file,
  genotypes = c("AA", "AB", "BB"),
  na.strings = c("NA"),
  estimate.map = FALSE
))

cross_phys <- ensure_map_names(cross_phys)
cross_phys$pheno <- data.frame(
  MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
  row.names = common_ids,
  stringsAsFactors = FALSE
)

phe_col <- "MeanDiseaseIncidence"

cross_phys <- run_step("Jittering map", {
  jittermap(cross_phys)
})
cross_phys <- ensure_map_names(cross_phys)

cross_phys <- run_step("Dropping null markers", {
  drop.nullmarkers(cross_phys)
})
cross_phys <- ensure_map_names(cross_phys)

safe_write_capture(
  summary(cross_phys),
  file = file.path(output_dir, "cross_phys_summary.txt")
)

# Keep the true F2 object for downstream QTL/MQM work.
cross_qtl_f2 <- ensure_map_names(cross_phys)
cross_qtl_f2$pheno <- data.frame(
  MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
  row.names = common_ids,
  stringsAsFactors = FALSE
)

# -----------------------------
# Early linkage diagnostics on the physical cross
# -----------------------------
cross_phys <- run_step("Estimating pairwise recombination fractions", {
  suppressWarnings(est.rf(cross_phys))
})

allele_swaps <- tryCatch(
  checkAlleles(cross_phys, threshold = 5),
  error = function(e) NULL
)

if (!is.null(allele_swaps) && nrow(allele_swaps) > 0) {
  write.csv(
    allele_swaps,
    file = file.path(output_dir, "CheckAlleles_Suspects.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  swap_markers <- if ("marker" %in% names(allele_swaps)) allele_swaps$marker else allele_swaps[[1]]
  swap_markers <- keep_nonempty_chr(swap_markers)
  
  if (length(swap_markers) > 0) {
    cross_phys <- run_step("Switching likely allele-flipped markers", {
      switchAlleles(cross_phys, swap_markers)
    })
    
    cross_phys <- run_step("Re-estimating recombination fractions after allele switch", {
      suppressWarnings(est.rf(cross_phys))
    })
    
    cross_qtl_f2 <- ensure_map_names(cross_phys)
    cross_qtl_f2$pheno <- data.frame(
      MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
      row.names = common_ids,
      stringsAsFactors = FALSE
    )
  }
}

# -----------------------------
# ASMap genetic map construction
# -----------------------------
cross_gen <- NULL
asmap_ok <- FALSE

if (use_asmap) {
  chr_levels <- sort_chr_labels(names(cross_phys$geno))
  
  cross_bcsft <- tryCatch(
    convert2bcsft(
      cross_phys,
      BC.gen = asmap_BC_gen,
      F.gen = asmap_F_gen,
      estimate.map = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(cross_bcsft)) {
    if (is.null(rownames(cross_bcsft$pheno))) {
      rownames(cross_bcsft$pheno) <- common_ids[seq_len(nrow(cross_bcsft$pheno))]
    }
    cross_bcsft$pheno$ID <- rownames(cross_bcsft$pheno)
    
    cross_gen <- tryCatch(
      run_step("ASMap MSTmap construction", {
        cross_tmp <- NULL
        invisible(capture.output(
          cross_tmp <- ASMap::mstmap(
            cross_bcsft,
            chr = chr_levels,
            id = "ID",
            bychr = TRUE,
            suffix = "numeric",
            anchor = FALSE,
            dist.fun = "kosambi",
            objective.fun = "COUNT",
            p.value = asmap_p_value,
            noMap.dist = asmap_noMap_dist,
            noMap.size = asmap_noMap_size,
            miss.thresh = asmap_miss_thresh,
            mvest.bc = FALSE,
            detectBadData = FALSE,
            return.imputed = FALSE,
            trace = asmap_trace
          )
        ))
        cross_tmp
      }),
      error = function(e) NULL
    )
    
    if (!is.null(cross_gen)) {
      cross_gen <- ensure_map_names(cross_gen)
      
      if (length(names(cross_gen$geno)) != 17 || any(grepl("\\.", names(cross_gen$geno)))) {
        cross_gen <- tryCatch(
          run_step("Merging ASMap subgroups back to chromosome groups", {
            merge_asmap_to_base_chromosomes(cross_gen)
          }),
          error = function(e) cross_gen
        )
        cross_gen <- ensure_map_names(cross_gen)
      }
      
      cross_gen$pheno <- data.frame(
        MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
        row.names = common_ids,
        stringsAsFactors = FALSE
      )
      
      cross_gen <- run_step("Jittering ASMap map", {
        jittermap(cross_gen)
      })
      cross_gen <- ensure_map_names(cross_gen)
      
      cross_gen <- run_step("Dropping null markers from ASMap map", {
        drop.nullmarkers(cross_gen)
      })
      cross_gen <- ensure_map_names(cross_gen)
      
      asmap_ok <- TRUE
      
      safe_write_capture(
        summary(cross_gen),
        file = file.path(output_dir, "cross_asmap_summary.txt")
      )
      
      safe_write_capture(
        summaryMap(cross_gen),
        file = file.path(output_dir, "ASMap_Genetic_Map_Summary.txt")
      )
    }
  }
}

if (!asmap_ok) {
  cross_gen <- cross_phys
  safe_write_capture(
    summaryMap(cross_gen),
    file = file.path(output_dir, "Physical_Fallback_Map_Summary.txt")
  )
}

# -----------------------------
# QTL/MQM analysis object
# Keep it F2. Apply the ASMap map only if it can be aligned cleanly.
# -----------------------------
cross_qtl_for_analysis <- ensure_map_names(cross_qtl_f2)
cross_qtl_for_analysis$pheno <- data.frame(
  MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
  row.names = common_ids,
  stringsAsFactors = FALSE
)

if (asmap_ok) {
  aligned_map <- tryCatch(
    align_asmap_map_to_f2(cross_qtl_f2, cross_gen),
    error = function(e) NULL
  )
  
  if (!is.null(aligned_map)) {
    cross_qtl_for_analysis <- apply_map_list_to_cross(cross_qtl_for_analysis, aligned_map)
    cross_qtl_for_analysis$pheno <- data.frame(
      MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
      row.names = common_ids,
      stringsAsFactors = FALSE
    )
  }
}

if (!is_f2_cross(cross_qtl_for_analysis)) {
  cross_qtl_for_analysis <- ensure_map_names(cross_qtl_f2)
  cross_qtl_for_analysis$pheno <- data.frame(
    MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
    row.names = common_ids,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Map tables
# -----------------------------
final_map_tbl <- build_map_table(cross_gen, marker_lookup, pos_label = "cM_final")

write.csv(
  final_map_tbl,
  file = file.path(output_dir, "ASMap_Genetic_Map_Table.csv"),
  row.names = FALSE,
  quote = FALSE
)

subgroup_summary <- final_map_tbl %>%
  count(chr_main, chr_label, name = "n_markers") %>%
  arrange(chr_main, chr_label)

write.csv(
  subgroup_summary,
  file = file.path(output_dir, "ASMap_Chromosome_Subgroup_Summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

base_chr_summary <- final_map_tbl %>%
  count(chr_main, name = "n_markers") %>%
  arrange(chr_main)

write.csv(
  base_chr_summary,
  file = file.path(output_dir, "ASMap_Chromosome_Summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# QC plots
# -----------------------------
cross_err <- run_step("Calculating genotyping error LODs", {
  calc.errorlod(cross_phys, error.prob = error_prob)
})

toperr <- tryCatch(
  top.errorlod(cross_err, cutoff = 6),
  error = function(e) NULL
)

if (!is.null(toperr) && nrow(toperr) > 0) {
  write.csv(
    toperr,
    file = file.path(output_dir, "Top_ErrorLOD_Genotypes.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

safe_write_capture(
  summary(cross_qtl_for_analysis),
  file = file.path(output_dir, "cross_qtl_summary.txt")
)

safe_pdf(file.path(output_dir, "R_qtl_QC_and_Diagnostics.pdf"), width = 11, height = 8.5, {
  par(mfrow = c(2, 2), las = 1)
  plotMissing(cross_phys, main = "Missing genotype pattern")
  plot(ntyped(cross_phys), ylab = "No. typed markers", main = "No. genotypes by individual")
  plot(ntyped(cross_phys, "mar"), ylab = "No. typed individuals", main = "No. genotypes by marker")
  
  gt <- geno.table(cross_phys, scanone.output = TRUE)
  plot(gt, ylab = expression(paste(-log[10], " P-value")), main = "Segregation distortion")
  
  par(mfrow = c(1, 2), las = 1)
  cg <- comparegeno(cross_phys)
  hist(
    cg[lower.tri(cg)],
    breaks = seq(0, 1, len = 101),
    xlab = "Proportion matching genotypes",
    main = "Pairwise genotype similarity"
  )
  plot(countXO(cross_phys), ylab = "Number of crossovers", main = "Crossovers per individual")
  
  par(mfrow = c(1, 1))
  plotMap(cross_phys, show.marker.names = FALSE, main = "Physical-position map (Mb)")
})

# -----------------------------
# QTL scans
# -----------------------------
cross_qtl_for_analysis <- run_step("Precomputing genotype probabilities", {
  calc.genoprob(
    cross_qtl_for_analysis,
    step = 1,
    error.prob = error_prob,
    map.function = "kosambi"
  )
})

scan_result_mr <- run_step("Marker regression scan", {
  scanone(cross_qtl_for_analysis, pheno.col = phe_col, method = "mr")
})

scan_result_kw <- run_step("Manual Kruskal-Wallis scan", {
  manual_kruskal_scan(cross_qtl_for_analysis, phe_col)
})

scan_result_hk <- run_step("Haley-Knott interval mapping scan", {
  scanone(cross_qtl_for_analysis, pheno.col = phe_col, method = "hk")
})

scan_result_em <- run_step("EM interval mapping scan", {
  scanone(cross_qtl_for_analysis, pheno.col = phe_col, method = "em")
})

perm_mr <- run_step("Permutation test for marker regression", {
  scanone(cross_qtl_for_analysis, pheno.col = phe_col, method = "mr", n.perm = perm_n, verbose = TRUE)
})

perm_em <- run_step("Permutation test for interval mapping", {
  scanone(cross_qtl_for_analysis, pheno.col = phe_col, method = "em", n.perm = perm_n, verbose = TRUE)
})

threshold_mr <- as.numeric(summary(perm_mr, alpha = 0.05))
threshold_em <- as.numeric(summary(perm_em, alpha = 0.05))

writeLines(
  paste0("Genome-wide 0.05 threshold (marker regression): ", threshold_mr),
  con = file.path(output_dir, "Permutation_Threshold_MR.txt")
)

writeLines(
  paste0("Genome-wide 0.05 threshold (EM interval mapping): ", threshold_em),
  con = file.path(output_dir, "Permutation_Threshold_EM.txt")
)

# -----------------------------
# QTL plots
# -----------------------------
png(file.path(output_dir, "SingleQTL_MarkerRegression.png"),
    width = 1800, height = 1200, res = 200)
plot(scan_result_mr, main = "Single-QTL Scan (Marker Regression)")
abline(h = threshold_mr, col = "red", lty = 2)
dev.off()

p_kw <- ggplot(
  scan_result_kw %>% filter(!is.na(pos), !is.na(kw_score)),
  aes(x = pos, y = kw_score)
) +
  geom_point(size = 0.7) +
  facet_wrap(~ chr, scales = "free_x", ncol = 4) +
  labs(
    title = "Manual Kruskal-Wallis scan",
    x = "Position",
    y = expression(-log[10](p))
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "SingleQTL_KruskalWallis.png"),
  plot = p_kw,
  width = 14,
  height = 10,
  dpi = 200
)

safe_pdf(file.path(output_dir, "SingleQTL_KruskalWallis.pdf"), width = 14, height = 10, {
  print(p_kw)
})

png(file.path(output_dir, "SingleQTL_HK.png"),
    width = 1800, height = 1200, res = 200)
plot(scan_result_hk, main = "Single-QTL Scan (Haley-Knott)")
dev.off()

png(file.path(output_dir, "SingleQTL_EM.png"),
    width = 1800, height = 1200, res = 200)
plot(scan_result_em, main = "Single-QTL Scan (EM interval mapping)")
abline(h = threshold_em, col = "red", lty = 2)
dev.off()

safe_pdf(file.path(output_dir, "R_qtl_QTL_Scans.pdf"), width = 11, height = 8.5, {
  plot(scan_result_mr, main = "Single-QTL Scan (Marker Regression)")
  abline(h = threshold_mr, col = "red", lty = 2)
  
  print(p_kw)
  
  plot(scan_result_hk, main = "Single-QTL Scan (Haley-Knott)")
  
  plot(scan_result_em, main = "Single-QTL Scan (EM Interval Mapping)")
  abline(h = threshold_em, col = "red", lty = 2)
  
  plot(perm_mr, main = "Permutation test: Marker Regression")
  plot(perm_em, main = "Permutation test: EM Interval Mapping")
})

# -----------------------------
# Primary QTL peak selection
# -----------------------------
primary_peak_mr <- as.data.frame(scan_result_mr[which.max(scan_result_mr$lod), , drop = FALSE])
primary_peak_mr$chr <- as.character(primary_peak_mr$chr)
primary_peak_mr$chr_main <- chr_main_from_label(primary_peak_mr$chr)
primary_peak_mr$pos <- as.numeric(primary_peak_mr$pos)
primary_peak_mr$lod <- as.numeric(primary_peak_mr$lod)
primary_peak_mr$peak_label <- paste0("Chr ", primary_peak_mr$chr_main, " @ ", round(primary_peak_mr$pos, 2), " cM")

write.csv(
  primary_peak_mr,
  file = file.path(output_dir, "Primary_MarkerRegression_Peak.csv"),
  row.names = FALSE,
  quote = FALSE
)

qtl_peaks_em <- tryCatch(
  summary(scan_result_em, perms = perm_em, alpha = 0.05),
  error = function(e) NULL
)

if (is.null(qtl_peaks_em) || nrow(as.data.frame(qtl_peaks_em)) == 0) {
  qtl_peaks_em <- as.data.frame(scan_result_em[which.max(scan_result_em$lod), , drop = FALSE])
} else {
  qtl_peaks_em <- as.data.frame(qtl_peaks_em)
}

qtl_peaks_em$chr <- as.character(qtl_peaks_em$chr)
qtl_peaks_em$chr_main <- chr_main_from_label(qtl_peaks_em$chr)
qtl_peaks_em$pos <- as.numeric(qtl_peaks_em$pos)
qtl_peaks_em$lod <- as.numeric(qtl_peaks_em$lod)
qtl_peaks_em$peak_label <- paste0("Chr ", qtl_peaks_em$chr_main, " @ ", round(qtl_peaks_em$pos, 2), " cM")

write.csv(
  qtl_peaks_em,
  file = file.path(output_dir, "Significant_QTL_Peaks_EM.csv"),
  row.names = FALSE,
  quote = FALSE
)

top_lod_sites <- as.data.frame(scan_result_em[order(-scan_result_em$lod), , drop = FALSE])
top_lod_sites <- head(top_lod_sites, 20)
top_lod_sites$chr <- as.character(top_lod_sites$chr)
top_lod_sites$chr_main <- chr_main_from_label(top_lod_sites$chr)

write.csv(
  top_lod_sites,
  file = file.path(output_dir, "Top_LOD_Sites_EM.csv"),
  row.names = FALSE,
  quote = FALSE
)

chr_lod_summary <- as.data.frame(scan_result_em) %>%
  mutate(chr_main = chr_main_from_label(as.character(chr))) %>%
  group_by(chr_main) %>%
  summarise(
    max_lod = max(lod, na.rm = TRUE),
    best_pos = pos[which.max(lod)][1],
    .groups = "drop"
  ) %>%
  arrange(desc(max_lod))

write.csv(
  chr_lod_summary,
  file = file.path(output_dir, "Chromosome_MaxLOD_Summary_EM.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  as.data.frame(scan_result_mr),
  file = file.path(output_dir, "ScanOne_MarkerRegression_Results.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  scan_result_kw,
  file = file.path(output_dir, "Manual_KruskalWallis_Results.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  as.data.frame(scan_result_hk),
  file = file.path(output_dir, "ScanOne_HaleyKnott_Results.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  as.data.frame(scan_result_em),
  file = file.path(output_dir, "ScanOne_EM_Results.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# QTL peak genotype-effect analysis
# -----------------------------
primary_peak_tbl <- primary_peak_mr

qtl_effects <- extract_qtl_peak_effects(
  cross_obj = cross_qtl_for_analysis,
  peaks_tbl = primary_peak_tbl,
  map_tbl = final_map_tbl,
  pheno_col = "MeanDiseaseIncidence"
)

if (nrow(qtl_effects$effect) > 0) {
  write.csv(
    qtl_effects$effect,
    file = file.path(output_dir, "QTL_Peak_Genotype_Effects.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

if (nrow(qtl_effects$quartile) > 0) {
  write.csv(
    qtl_effects$quartile,
    file = file.path(output_dir, "QTL_Peak_Quartile_Association.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

if (nrow(qtl_effects$plot_df) > 0) {
  p_qtl_assoc <- ggplot(qtl_effects$plot_df, aes(x = genotype_code, y = phenotype)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    facet_wrap(~ peak_label, scales = "free_x") +
    labs(
      title = "Phenotype association at the primary QTL peak",
      x = "Genotype code",
      y = "Mean disease incidence"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    file.path(output_dir, "QTL_Peak_Genotype_Association.png"),
    plot = p_qtl_assoc,
    width = 10,
    height = 7,
    dpi = 200
  )
  
  safe_pdf(file.path(output_dir, "QTL_Peak_Genotype_Association.pdf"), width = 10, height = 7, {
    print(p_qtl_assoc)
  })
}

# -----------------------------
# MQM mapping
# -----------------------------
cross_mqm_source <- cross_qtl_for_analysis
if (!is_f2_cross(cross_mqm_source)) {
  cross_mqm_source <- ensure_map_names(cross_qtl_f2)
  cross_mqm_source$pheno <- data.frame(
    MeanDiseaseIncidence = as.numeric(pheno_filt$MeanDiseaseIncidence),
    row.names = common_ids,
    stringsAsFactors = FALSE
  )
}

mqm_out <- run_mqm_analysis(cross_mqm_source, phe_col)

if (!is.null(mqm_out)) {
  saveRDS(mqm_out$cofactors, file = file.path(output_dir, "mqm_cofactors.rds"))
  saveRDS(mqm_out$result, file = file.path(output_dir, "mqm_result.rds"))
  
  if (inherits(mqm_out$result, "scanone")) {
    mqm_df <- as.data.frame(mqm_out$result)
    write.csv(
      mqm_df,
      file = file.path(output_dir, "MQM_Scan_Results.csv"),
      row.names = FALSE,
      quote = FALSE
    )
    
    png(file.path(output_dir, "MQM_Scan.png"),
        width = 1800, height = 1200, res = 200)
    plot(mqm_out$result, main = "MQM scan")
    dev.off()
  }
  
  writeLines(
    paste0("MQM completed using ", mqm_out$ncof, " cofactors."),
    con = file.path(output_dir, "MQM_Status.txt")
  )
} else {
  writeLines(
    "MQM was skipped because no stable model could be fit without errors.",
    con = file.path(output_dir, "MQM_Status.txt")
  )
}

# -----------------------------
# Marker annotation tables
# -----------------------------
marker_annot <- parse_marker_annotations(all_markers_file, marker_9k_file)

marker_annot <- marker_annot %>%
  left_join(marker_lookup, by = c("join_name" = "SNPName"))

marker_annot <- marker_annot %>%
  left_join(
    final_map_tbl %>% select(marker_id, map_chr_main = chr_main, cM_final),
    by = "marker_id"
  ) %>%
  mutate(
    chr_main = coalesce(
      map_chr_main,
      chr_final,
      chr_num_9k,
      suppressWarnings(as.integer(stringr::str_extract(linkage_group, "\\d+")))
    )
  )

write.csv(
  marker_annot,
  file = file.path(output_dir, "Marker_Annotation_Merged.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Candidate markers from the primary QTL peak
# -----------------------------
marker_annot_mapped <- marker_annot %>%
  filter(!is.na(marker_id), !is.na(chr_main), !is.na(cM_final))

peak_chr_main <- primary_peak_tbl$chr_main
peak_pos <- primary_peak_tbl$pos

cand_flag <- vapply(seq_len(nrow(marker_annot_mapped)), function(i) {
  chr_i <- marker_annot_mapped$chr_main[i]
  cm_i <- marker_annot_mapped$cM_final[i]
  if (is.na(chr_i) || is.na(cm_i)) return(FALSE)
  
  peak_pos_chr <- peak_pos[peak_chr_main == chr_i]
  peak_pos_chr <- peak_pos_chr[!is.na(peak_pos_chr)]
  if (length(peak_pos_chr) == 0) return(FALSE)
  
  any(abs(cm_i - peak_pos_chr) <= qtl_window_cm)
}, logical(1))

candidate_markers <- marker_annot_mapped[cand_flag, , drop = FALSE] %>%
  arrange(chr_main, cM_final, marker_id)

if (nrow(candidate_markers) == 0) {
  top_peak <- primary_peak_tbl[1, , drop = FALSE]
  candidate_markers <- marker_annot_mapped %>%
    filter(chr_main == top_peak$chr_main) %>%
    mutate(dist_to_peak = abs(cM_final - top_peak$pos)) %>%
    arrange(dist_to_peak, chr_main, cM_final, marker_id) %>%
    slice_head(n = candidate_keep_top_n_if_none)
  
  if (nrow(candidate_markers) == 0) {
    candidate_markers <- marker_annot_mapped %>%
      mutate(dist_to_peak = abs(cM_final - top_peak$pos)) %>%
      arrange(dist_to_peak, chr_main, cM_final, marker_id) %>%
      slice_head(n = candidate_keep_top_n_if_none)
  }
}

write.csv(
  candidate_markers,
  file = file.path(output_dir, "Candidate_Markers_For_Figures.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Save workspace objects
# -----------------------------
saveRDS(cross_phys, file = file.path(output_dir, "cross_phys_final.rds"))
saveRDS(cross_qtl_f2, file = file.path(output_dir, "cross_qtl_f2_final.rds"))
saveRDS(cross_qtl_for_analysis, file = file.path(output_dir, "cross_qtl_for_analysis_final.rds"))
saveRDS(scan_result_mr, file = file.path(output_dir, "scan_result_mr.rds"))
saveRDS(scan_result_kw, file = file.path(output_dir, "scan_result_kw.rds"))
saveRDS(scan_result_hk, file = file.path(output_dir, "scan_result_hk.rds"))
saveRDS(scan_result_em, file = file.path(output_dir, "scan_result_em.rds"))
saveRDS(perm_mr, file = file.path(output_dir, "perm_scan_mr.rds"))
saveRDS(perm_em, file = file.path(output_dir, "perm_scan_em.rds"))
saveRDS(primary_peak_mr, file = file.path(output_dir, "primary_peak_mr.rds"))
saveRDS(qtl_peaks_em, file = file.path(output_dir, "qtl_peaks_em.rds"))
saveRDS(marker_annot, file = file.path(output_dir, "marker_annot_merged.rds"))
saveRDS(physical_map_tbl, file = file.path(output_dir, "physical_map_tbl.rds"))
saveRDS(final_map_tbl, file = file.path(output_dir, "genetic_map_tbl.rds"))
saveRDS(qtl_effects, file = file.path(output_dir, "qtl_effects.rds"))

message("Done. Outputs are in: ", output_dir)