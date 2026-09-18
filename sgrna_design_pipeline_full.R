############################################################
# Configurable sgRNA design and interactive report pipeline
#
# Supports crisprDesign guide generation or existing RDS input.
# Edit default_config() before running, or pass an .R/.rds config to Rscript.
# Required choices are grouped at the top of the file:
# genome = BSgenome package, FASTA, or raw sequence
# target = gene (Entrez ID), genomic region, or custom sequence
# annotation = TxDb package/object and transcript selection
# nuclease = built-in/custom PAM and spacer definition
# modality = CRISPRko/i/a, base_edit, or prime_edit
# filters = score, GC, homopolymer, strand, and custom filters
# off_target = optional Bowtie alignment/scoring
# paired_guides/editing/cloning = optional advanced design modules
# output = CSV/RDS/HTML/BED/pairs/transcript-map destinations
# The report uses 1-based mature-transcript (spliced cDNA) positions.
############################################################

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

default_config <- function() {
 list(
 mode = "design", # design or rds
 project_name = "TP53_sgRNA_design",
 runtime = list(install_missing = FALSE, quiet = FALSE),
 genome = list(
 source = "bsgenome", # bsgenome, fasta, or sequence
 package = "BSgenome.Hsapiens.UCSC.hg38",
 object = "BSgenome.Hsapiens.UCSC.hg38",
 fasta = NULL,
 sequence = NULL
 ),
 target = list(
 type = "gene", # region, gene, or sequence
 name = "TP53",
 gene_id = "7157",
 gene_column = "gene_id",
 feature_type = "cds",
 seqnames = NULL,
 start = NULL,
 end = NULL,
 sequence = NULL,
 sequence_start = NULL
 ),
 annotation = list(
 enabled = TRUE,
 txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene",
 txdb_object = "TxDb.Hsapiens.UCSC.hg38.knownGene",
 organism_package = "org.Hs.eg.db",
 organism_keytype = "ENTREZID"
 ),
 # Transcript-relative coordinates are reported on the mature spliced
 # transcript (1-based cDNA position). Leave selected_id = NULL to
 # show the first annotated transcript, or set a tx_name explicitly.
 transcript = list(
 coordinate_system = "mature_transcript",
 selected_id = NULL,
 selector_in_report = TRUE
 ),
 # CRISPRko, CRISPRi, CRISPRa, base_edit, or prime_edit. CRISPRi/a
 # optionally restricts guides to a strand-aware TSS window.
 modality = list(
 type = "CRISPRko",
 tss_window = c(NA_real_, NA_real_)
 ),
 cell_context = list(
 name = "unspecified", species = NULL, strain = NULL,
 vcf = NULL, metadata = list()
 ),
 nuclease = list(
 name = "SpCas9", pams = NULL, weights = NULL,
 pam_side = "3prime", spacer_length = 20L, target_type = "DNA",
 metadata = list()
 ),
 design = list(
 canonical = TRUE, both_strands = TRUE, strict_overlap = TRUE,
 remove_ambiguities = TRUE, remove_duplicates = TRUE
 ),
 scoring = list(
 # These two methods are fast, reproducible, and available for SpCas9.
 # Set add_on_target = FALSE only when supplying pre-scored guides.
 add_on_target = TRUE,
 methods = c("ruleset1", "crisprater"),
 add_sequence_features = TRUE, add_composite = TRUE
 ),
 off_target = list(
 enabled = FALSE, bowtie_index = NULL, mismatches = 3L,
 all_alignments = TRUE
 ),
 paired_guides = list(
 enabled = TRUE,
 min_distance = 20L,
 max_distance = 1000L,
 opposite_strands_only = FALSE
 ),
 editing = list(
 # For base editing, positions are 1-based in the protospacer.
 window = 4:8,
 from = "C",
 to = "T",
 pbs_length = 13L,
 rtt_length = 20L,
 pbs_sequence = NULL,
 rtt_sequence = NULL
 ),
 cloning = list(
 enabled = FALSE,
 five_prime_adapter = "",
 three_prime_adapter = "",
 scaffold = "",
 restriction_sites = character()
 ),
 filters = list(
 min_composite_score = NA_real_, min_gc = NA_real_, max_gc = NA_real_,
 max_homopolymer = 4L, strands = c("+", "-"),
 exclude_sequences = character(), custom = NULL
 ),
 input = list(rds_file = NULL),
 output = list(
 directory = "/home/goro/R/output", prefix = "TP53_sgRNA_design",
 write_rds = TRUE, write_csv = TRUE, write_html = TRUE,
 write_bed = TRUE, write_pairs = TRUE, write_transcript_map = TRUE,
 top_guides = 20L
 )
 )
}

modify_config <- function(config, updates) {
 if (!is.list(config) || !is.list(updates)) stop("config and updates must be lists")
 merge_lists <- function(x, y) {
 for (nm in names(y)) {
 if (is.list(y[[nm]]) && is.list(x[[nm]] %||% NULL)) x[[nm]] <- merge_lists(x[[nm]], y[[nm]])
 else x[[nm]] <- y[[nm]]
 }
 x
 }
 merge_lists(config, updates)
}

load_config_file <- function(path, base = default_config()) {
 if (!file.exists(path)) stop("Configuration file not found: ", path)
 user_config <- if (tolower(tools::file_ext(path)) == "rds") {
 readRDS(path)
 } else {
 env <- new.env(parent = globalenv())
 sys.source(path, envir = env)
 env$config %||% env$CONFIG
 }
 if (!is.list(user_config)) stop("Config file must define a list named 'config' or 'CONFIG'")
 modify_config(base, user_config)
}

validate_config <- function(config) {
 needed <- c("mode", "genome", "target", "nuclease", "design", "filters", "output")
 missing <- setdiff(needed, names(config))
 if (length(missing)) stop("Missing configuration sections: ", paste(missing, collapse = ", "))
 config$mode <- match.arg(config$mode, c("design", "rds"))
 config$genome$source <- match.arg(config$genome$source, c("bsgenome", "fasta", "sequence"))
 config$target$type <- match.arg(config$target$type, c("region", "gene", "sequence"))
 config$transcript$coordinate_system <- match.arg(config$transcript$coordinate_system, "mature_transcript")
 config$modality$type <- match.arg(config$modality$type, c("CRISPRko", "CRISPRi", "CRISPRa", "base_edit", "prime_edit"))
 if (config$mode == "rds" && is.null(config$input$rds_file)) stop("input$rds_file is required for mode = 'rds'")
 if (config$target$type == "region") {
 if (is.null(config$target$seqnames) || is.null(config$target$start) || is.null(config$target$end))
 stop("A region target requires seqnames, start, and end")
 if (config$target$start > config$target$end) stop("target start must be <= target end")
 }
 if (config$target$type == "gene" && config$genome$source != "bsgenome") {
 stop("gene targets require a BSgenome source so genomic coordinates can be resolved")
 }
 if (config$genome$source == "fasta" && is.null(config$genome$fasta)) stop("genome$fasta is required for FASTA input")
 if (config$genome$source == "sequence" && is.null(config$genome$sequence) && is.null(config$target$sequence))
 stop("Provide genome$sequence or target$sequence for sequence input")
 if (is.null(config$output$directory) || is.null(config$output$prefix)) stop("output directory and prefix are required")
 if (length(config$modality$tss_window) != 2L) stop("modality$tss_window must contain two values")
 if (isTRUE(config$paired_guides$enabled) && config$paired_guides$min_distance > config$paired_guides$max_distance)
 stop("paired_guides$min_distance must be <= max_distance")
 if (config$modality$type == "base_edit") {
 if (any(config$editing$window < 1L)) stop("editing$window positions must be positive")
 if (nchar(config$editing$from) != 1L || nchar(config$editing$to) != 1L)
 stop("editing$from and editing$to must be single bases")
 }
 config
}

install_required_packages <- function(config) {
 required <- c("crisprBase", "crisprDesign", "jsonlite")
 optional <- character()
 if (config$genome$source == "bsgenome") optional <- c(optional, config$genome$package)
 if (config$genome$source == "fasta") optional <- c(optional, "Rsamtools", "Biostrings")
 if (isTRUE(config$annotation$enabled)) optional <- c(optional, "GenomicFeatures", config$annotation$txdb_package)
 if (isTRUE(config$annotation$enabled) && !is.null(config$annotation$organism_package)) optional <- c(optional, config$annotation$organism_package)
 all_packages <- unique(c(required, optional))
 missing <- all_packages[!vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)]
 if (!length(missing)) return(invisible(TRUE))
 if (!isTRUE(config$runtime$install_missing)) stop("Missing packages: ", paste(missing, collapse = ", "), ". Install them or set runtime$install_missing = TRUE.")
 if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
 bioc <- intersect(missing, c("crisprBase", "crisprDesign", "GenomicFeatures", "Rsamtools", "Biostrings", optional))
 if (length(bioc)) BiocManager::install(bioc, ask = FALSE, update = FALSE)
 cran <- setdiff(missing, bioc)
 if (length(cran)) install.packages(cran, repos = "https://cloud.r-project.org")
 invisible(TRUE)
}

resolve_nuclease <- function(config) {
 name <- config$nuclease$name
 if (is.null(config$nuclease$pams)) {
 available <- crisprBase::getAvailableCrisprNucleases()
 if (!(name %in% available)) stop("Unknown nuclease '", name, "'. Available: ", paste(available, collapse = ", "), ". Provide nuclease$pams for custom motifs.")
 return(getExportedValue("crisprBase", name))
 }
 do.call(crisprBase::CrisprNuclease, list(
 nucleaseName = name, pams = config$nuclease$pams, weights = config$nuclease$weights,
 pam_side = config$nuclease$pam_side, spacer_length = config$nuclease$spacer_length,
 targetType = config$nuclease$target_type, metadata = config$nuclease$metadata
 ))
}

resolve_bsgenome <- function(config) {
 if (config$genome$source != "bsgenome") return(NULL)
 if (is.null(config$genome$object)) stop("genome$object must name an exported BSgenome object")
 if (!requireNamespace(config$genome$package, quietly = TRUE)) stop("Genome package is not installed: ", config$genome$package)
 getExportedValue(config$genome$package, config$genome$object)
}

resolve_txdb <- function(config) {
 if (!isTRUE(config$annotation$enabled) || is.null(config$annotation$txdb_package)) return(NULL)
 if (!requireNamespace(config$annotation$txdb_package, quietly = TRUE)) return(NULL)
 # Some Bioconductor installations keep the TxDb SQLite file in a
 # read-only package directory. Gene queries may need a temporary lock,
 # so load a writable copy when the packaged object cannot be queried.
 if (identical(config$target$type, "gene") && requireNamespace("GenomicFeatures", quietly = TRUE)) {
 db_file <- system.file(
 "extdata", paste0(config$annotation$txdb_object, ".sqlite"),
 package = config$annotation$txdb_package
 )
 if (nzchar(db_file) && file.exists(db_file)) {
 tmp_file <- file.path(tempdir(), paste0(config$annotation$txdb_object, "_", Sys.getpid(), ".sqlite"))
 if (!file.exists(tmp_file)) file.copy(db_file, tmp_file, overwrite = TRUE)
 return(AnnotationDbi::loadDb(tmp_file))
 }
 }
 getExportedValue(config$annotation$txdb_package, config$annotation$txdb_object)
}

resolve_sequence <- function(config, bsgenome) {
 if (!is.null(config$target$sequence)) return(as.character(config$target$sequence)[1L])
 if (config$genome$source == "sequence") return(as.character(config$genome$sequence)[1L])
 if (config$genome$source == "bsgenome") return(as.character(Biostrings::getSeq(bsgenome, names = config$target$seqnames, start = config$target$start, end = config$target$end)))
 fa <- Rsamtools::FaFile(config$genome$fasta)
 gr <- GenomicRanges::GRanges(config$target$seqnames, IRanges::IRanges(config$target$start, config$target$end))
 as.character(Rsamtools::scanFa(fa, gr)[[1L]])
}

target_granges <- function(config) {
 GenomicRanges::GRanges(config$target$seqnames, IRanges::IRanges(config$target$start, config$target$end), strand = "*")
}

first_column <- function(data, candidates) {
 hit <- candidates[candidates %in% colnames(data)]
 if (length(hit)) data[[hit[1L]]] else NULL
}

normalise_guides <- function(guides, config) {
 # GuideSet has a useful as.data.frame method. Older versions of
 # crisprDesign expose flattenGuideSet(), but it emits a deprecation
 # message and can return a diagnostic object instead of guide rows.
 table <- tryCatch(as.data.frame(guides, stringsAsFactors = FALSE), error = function(e) NULL)
 if (is.null(table) || ncol(table) <= 1L) {
 flat <- tryCatch(crisprDesign::flattenGuideSet(guides), error = function(e) guides)
 table <- as.data.frame(flat, stringsAsFactors = FALSE)
 }
 if (!nrow(table)) return(table)
 table$guide_id <- as.character(first_column(table, c("guide_id", "id", "ids")) %||% rownames(table))
 table$protospacer <- as.character(first_column(table, c("protospacer", "protospacers", "spacer", "sequence")) %||% NA_character_)
 table$pam <- as.character(first_column(table, c("pam", "pams", "PAM")) %||% NA_character_)
 table$strand <- as.character(first_column(table, c("strand", "spacer_strand")) %||% "*")
 table$seqnames <- as.character(first_column(table, c("seqnames", "seqname", "chromosome")) %||% config$target$seqnames %||% NA_character_)
 starts_raw <- first_column(table, c("start", "spacer_start", "protospacer_start"))
 ends_raw <- first_column(table, c("end", "spacer_end", "protospacer_end"))
 starts <- if (is.null(starts_raw)) rep(NA_integer_, nrow(table)) else suppressWarnings(as.integer(starts_raw))
 ends <- if (is.null(ends_raw)) rep(NA_integer_, nrow(table)) else suppressWarnings(as.integer(ends_raw))
 if ((!any(is.finite(starts)) || !any(is.finite(ends))) && !is.null(config$target$start)) {
 starts <- config$target$start + seq_len(nrow(table)) - 1L
 ends <- starts + config$nuclease$spacer_length - 1L
 }
 if (!is.null(config$target$sequence_start) && any(is.finite(starts)) && (config$genome$source == "sequence" || !is.null(config$target$sequence))) {
 starts <- starts + as.integer(config$target$sequence_start) - 1L
 ends <- ends + as.integer(config$target$sequence_start) - 1L
 }
 table$start <- starts; table$end <- ends
 cut_raw <- first_column(table, c("cut_site", "cutSite", "cut_sites"))
 cut <- if (is.null(cut_raw)) rep(NA_real_, nrow(table)) else suppressWarnings(as.numeric(cut_raw))
 table$cut_site <- cut
 if (all(is.na(cut)) && any(is.finite(starts))) table$cut_site <- starts
 score_columns <- grep("(^score|score$|on.?target|efficiency)", colnames(table), ignore.case = TRUE, value = TRUE)
 score_columns <- score_columns[vapply(table[score_columns], is.numeric, logical(1))]
 table$composite_score <- if (length(score_columns)) rowMeans(table[score_columns], na.rm = TRUE) else NA_real_
 table$gc_percent <- vapply(table$protospacer, function(s) {
 if (is.na(s) || !nzchar(s)) return(NA_real_)
 chars <- strsplit(toupper(s), "", fixed = TRUE)[[1L]]
 100 * sum(chars %in% c("G", "C")) / length(chars)
 }, numeric(1))
 table
}

try_add <- function(object, fun, ..., label) {
 tryCatch(fun(object, ...), error = function(e) { warning(label, " skipped: ", conditionMessage(e), call. = FALSE); object })
}

design_guides <- function(config, nuclease, bsgenome) {
 if (config$mode == "rds") {
 if (!file.exists(config$input$rds_file)) stop("RDS file not found: ", config$input$rds_file)
 return(readRDS(config$input$rds_file))
 }
 if (config$target$type == "gene") {
 txdb <- resolve_txdb(config)
 if (is.null(txdb)) stop("A gene target requires a usable TxDb annotation")
 # Query the TxDb directly. queryTxObject() expects a GRangesList and
 # can trigger a read-only lock when handed a TxDb object.
 feature_type <- config$target$feature_type
 by_gene <- switch(
 feature_type,
 cds = GenomicFeatures::cdsBy(txdb, by = "gene"),
 transcripts = GenomicFeatures::transcriptsBy(txdb, by = "gene"),
 exons = GenomicFeatures::exonsBy(txdb, by = "gene"),
 stop("Unsupported gene feature_type: ", feature_type,
 ". Use cds, transcripts, or exons.")
 )
 region <- by_gene[[as.character(config$target$gene_id)]]
 if (!length(region)) stop("No target features found for gene ", config$target$gene_id)
 return(crisprDesign::findSpacers(region, bsgenome = bsgenome, crisprNuclease = nuclease,
 canonical = config$design$canonical, both_strands = config$design$both_strands,
 strict_overlap = config$design$strict_overlap, remove_ambiguities = config$design$remove_ambiguities,
 remove_duplicates = config$design$remove_duplicates))
 }
 if (config$target$type == "sequence" || config$genome$source != "bsgenome") {
 return(crisprDesign::findSpacers(resolve_sequence(config, bsgenome), crisprNuclease = nuclease,
 both_strands = config$design$both_strands, remove_ambiguities = config$design$remove_ambiguities,
 remove_duplicates = config$design$remove_duplicates))
 }
 crisprDesign::findSpacers(target_granges(config), bsgenome = bsgenome, crisprNuclease = nuclease,
 canonical = config$design$canonical, both_strands = config$design$both_strands,
 strict_overlap = config$design$strict_overlap, remove_ambiguities = config$design$remove_ambiguities,
 remove_duplicates = config$design$remove_duplicates)
}

filter_guides <- function(table, config) {
 if (!nrow(table)) return(table)
 f <- config$filters; keep <- rep(TRUE, nrow(table))
 if (length(f$strands)) keep <- keep & (table$strand %in% f$strands | table$strand == "*")
 if (is.finite(f$min_composite_score)) keep <- keep & table$composite_score >= f$min_composite_score
 if (is.finite(f$min_gc)) keep <- keep & table$gc_percent >= f$min_gc
 if (is.finite(f$max_gc)) keep <- keep & table$gc_percent <= f$max_gc
 if (length(f$exclude_sequences)) keep <- keep & !(toupper(table$protospacer) %in% toupper(f$exclude_sequences))
 if (is.finite(f$max_homopolymer) && f$max_homopolymer > 1L) {
 pattern <- paste0("([ACGTN])\\1{", as.integer(f$max_homopolymer) - 1L, ",}")
 keep <- keep & !grepl(pattern, toupper(table$protospacer), perl = TRUE)
 }
 if (config$modality$type %in% c("CRISPRi", "CRISPRa") && all(is.finite(config$modality$tss_window)) && "distance_to_tss" %in% names(table)) {
 keep <- keep & table$distance_to_tss >= config$modality$tss_window[1L] & table$distance_to_tss <= config$modality$tss_window[2L]
 }
 if (is.function(f$custom)) {
 custom_keep <- f$custom(table)
 if (!is.logical(custom_keep) || length(custom_keep) != nrow(table)) stop("filters$custom must return one logical value per guide")
 keep <- keep & !is.na(custom_keep) & custom_keep
 }
 table[keep, , drop = FALSE]
}

get_gene_model <- function(config) {
 txdb <- resolve_txdb(config)
 if (is.null(txdb) || !requireNamespace("GenomicFeatures", quietly = TRUE) || !requireNamespace("GenomicRanges", quietly = TRUE)) return(NULL)
 tx <- if (identical(config$target$type, "gene")) {
 GenomicFeatures::transcriptsBy(txdb, by = "gene")[[as.character(config$target$gene_id)]]
 } else {
 GenomicFeatures::transcripts(txdb)
 }
 if (!is.null(config$target$seqnames) && !is.null(config$target$start)) {
 tx <- tx[as.character(GenomicRanges::seqnames(tx)) == config$target$seqnames & GenomicRanges::start(tx) <= config$target$end & GenomicRanges::end(tx) >= config$target$start]
 }
 if (!length(tx)) return(NULL)
 tx_table <- as.data.frame(tx, stringsAsFactors = FALSE); meta <- S4Vectors::mcols(tx)
 tx_ids <- as.character(meta$tx_id %||% seq_len(nrow(tx_table))); tx_table$tx_id <- tx_ids
 tx_table$transcript <- if ("tx_name" %in% colnames(meta)) as.character(meta$tx_name) else tx_ids
 tx_table$transcript <- make.unique(tx_table$transcript)
 tx_table$tss <- ifelse(as.character(tx_table$strand) == "+", tx_table$start, tx_table$end)
 exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = FALSE)
 exon_tables <- lapply(seq_len(nrow(tx_table)), function(i) {
 e <- exons[[tx_table$tx_id[i]]]; if (is.null(e) || !length(e)) return(NULL)
 out <- as.data.frame(e, stringsAsFactors = FALSE); out$tx_id <- tx_table$tx_id[i]; out$transcript <- tx_table$transcript[i]
 em <- S4Vectors::mcols(e); out$exon_number <- if ("exon_rank" %in% colnames(em)) as.integer(em$exon_rank) else seq_len(nrow(out)); out
 })
 list(transcripts = tx_table, exons = do.call(rbind, exon_tables), information = data.frame(SYMBOL = config$target$name, GENENAME = "", stringsAsFactors = FALSE))
}

empty_transcript_map <- function() {
 data.frame(
 guide_id = character(), transcript = character(), tx_id = character(),
 transcript_position = integer(), transcript_region = character(),
 exon_number = integer(), distance_to_tss = numeric(), strand = character(),
 composite_score = numeric(), seqnames = character(), cut_site = numeric(),
 stringsAsFactors = FALSE
 )
}

# Map genomic cleavage coordinates onto the mature, spliced transcript.
# Intronic cuts intentionally have NA transcript_position because introns
# are not part of the mature transcript; transcript_region records them.
annotate_transcript_positions <- function(table, model, config) {
 if (is.null(model) || !nrow(table) || !nrow(model$transcripts)) {
 return(list(table = table, map = empty_transcript_map()))
 }
 map_rows <- vector("list", nrow(table) * nrow(model$transcripts))
 k <- 0L
 for (tx_i in seq_len(nrow(model$transcripts))) {
 tx <- model$transcripts[tx_i, , drop = FALSE]
 ex <- model$exons[model$exons$tx_id == tx$tx_id, , drop = FALSE]
 if (!nrow(ex)) next
 ex <- ex[order(ex$start, decreasing = as.character(tx$strand) == "-"), , drop = FALSE]
 ex$width <- ex$end - ex$start + 1L
 ex$cumulative_start <- c(1L, head(cumsum(ex$width), -1L) + 1L)
 for (i in seq_len(nrow(table))) {
 cut <- suppressWarnings(as.numeric(table$cut_site[i]))
 same_chr <- !is.na(table$seqnames[i]) && as.character(table$seqnames[i]) == as.character(tx$seqnames)
 in_exon <- same_chr & cut >= ex$start & cut <= ex$end
 in_tx <- same_chr & cut >= min(tx$start, tx$end) & cut <= max(tx$start, tx$end)
 exon_number <- NA_integer_; tx_pos <- NA_integer_; region <- "outside_transcript"
 if (any(in_exon)) {
 j <- which(in_exon)[1L]
 tx_pos <- if (as.character(tx$strand) == "+") ex$cumulative_start[j] + cut - ex$start else ex$cumulative_start[j] + ex$end[j] - cut
 exon_number <- suppressWarnings(as.integer(ex$exon_number[j]))
 region <- "exon"
 } else if (in_tx) {
 region <- "intron"
 }
 tss_distance <- if (as.character(tx$strand) == "+") cut - tx$tss else tx$tss - cut
 k <- k + 1L
 map_rows[[k]] <- data.frame(
 guide_id = as.character(table$guide_id[i]), transcript = as.character(tx$transcript), tx_id = as.character(tx$tx_id),
 transcript_position = tx_pos, transcript_region = region, exon_number = exon_number,
 distance_to_tss = tss_distance, strand = as.character(table$strand[i]),
 composite_score = as.numeric(table$composite_score[i]), seqnames = as.character(table$seqnames[i]),
 cut_site = as.numeric(table$cut_site[i]), stringsAsFactors = FALSE
 )
 }
 }
 mapping <- if (k) do.call(rbind, map_rows[seq_len(k)]) else empty_transcript_map()
 selected <- config$transcript$selected_id %||% if (nrow(model$transcripts)) model$transcripts$transcript[1L] else NA_character_
 selected_rows <- mapping[mapping$transcript == selected, , drop = FALSE]
 if (!nrow(selected_rows) && nrow(mapping)) {
 selected <- mapping$transcript[1L]
 selected_rows <- mapping[mapping$transcript == selected, , drop = FALSE]
 }
 idx <- match(table$guide_id, selected_rows$guide_id)
 table$transcript <- selected
 table$mature_transcript_position <- selected_rows$transcript_position[idx]
 table$transcript_region <- selected_rows$transcript_region[idx]
 table$exon_number <- selected_rows$exon_number[idx]
 table$distance_to_tss <- selected_rows$distance_to_tss[idx]
 list(table = table, map = mapping)
}

annotate_modality <- function(table, config) {
 if (!nrow(table)) return(table)
 table$editing_candidates <- NA_character_
 if (config$modality$type == "base_edit") {
 w <- as.integer(config$editing$window); from <- toupper(config$editing$from); to <- toupper(config$editing$to)
 table$editing_candidates <- vapply(table$protospacer, function(seq) {
 chars <- strsplit(toupper(seq), "", fixed = TRUE)[[1L]]
 pos <- w[w <= length(chars) & chars[w] == from]
 if (!length(pos)) "" else paste0(from, pos, ">", to, collapse = ";")
 }, character(1))
 } else if (config$modality$type == "prime_edit") {
 table$editing_candidates <- paste0("PBS=", config$editing$pbs_length, ";RTT=", config$editing$rtt_length)
 }
 table
}

annotate_cloning <- function(table, config) {
 if (!nrow(table)) return(table)
 table$cloning_oligo <- NA_character_; table$restriction_site_hit <- NA_character_
 if (!isTRUE(config$cloning$enabled)) return(table)
 table$cloning_oligo <- paste0(config$cloning$five_prime_adapter, table$protospacer,
 config$cloning$three_prime_adapter, config$cloning$scaffold)
 if (length(config$cloning$restriction_sites)) {
 table$restriction_site_hit <- vapply(table$cloning_oligo, function(s) {
 hits <- config$cloning$restriction_sites[vapply(config$cloning$restriction_sites, function(m) grepl(m, s, fixed = TRUE), logical(1))]
 paste(hits, collapse = ";")
 }, character(1))
 }
 table
}

build_guide_pairs <- function(table, config) {
 if (!isTRUE(config$paired_guides$enabled) || nrow(table) < 2L) return(data.frame())
 rows <- list(); k <- 0L
 for (i in seq_len(nrow(table) - 1L)) for (j in (i + 1L):nrow(table)) {
 d <- abs(table$cut_site[j] - table$cut_site[i])
 if (!is.finite(d) || d < config$paired_guides$min_distance || d > config$paired_guides$max_distance) next
 if (isTRUE(config$paired_guides$opposite_strands_only) && table$strand[i] == table$strand[j]) next
 k <- k + 1L
 rows[[k]] <- data.frame(pair_id = paste0("pair_", k), guide_1 = table$guide_id[i], guide_2 = table$guide_id[j],
 cut_1 = table$cut_site[i], cut_2 = table$cut_site[j], deletion_bp = d,
 mean_score = mean(c(table$composite_score[i], table$composite_score[j]), na.rm = TRUE), stringsAsFactors = FALSE)
 }
 if (!k) return(data.frame())
 out <- do.call(rbind, rows); out[order(out$mean_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
}

write_bed <- function(table, path) {
 if (!nrow(table)) { file.create(path); return(invisible(path)) }
 bed <- data.frame(chrom = table$seqnames, start = pmax(0L, as.integer(table$cut_site) - 1L), end = as.integer(table$cut_site),
 name = table$guide_id, score = pmin(1000, pmax(0, round(1000 * table$composite_score))), strand = table$strand)
 utils::write.table(bed, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
 invisible(path)
}

escape_html <- function(x) {
 x <- gsub("&", "&amp;", as.character(x), fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE)
}

build_transcript_traces <- function(model) {
 if (is.null(model) || !nrow(model$transcripts)) return(list(traces = list(), names = character()))
 traces <- list()
 for (i in seq_len(nrow(model$transcripts))) {
 tx <- model$transcripts[i, , drop = FALSE]
 traces[[length(traces) + 1L]] <- list(type = "scatter", mode = "lines", x = c(tx$start, tx$end), y = c(i, i), line = list(color = "#7f8c8d", width = 2), showlegend = FALSE, hoverinfo = "text", text = paste0("Transcript: ", tx$transcript, "<br>", tx$seqnames, ":", tx$start, "-", tx$end, "<br>Strand: ", tx$strand))
 e <- model$exons[model$exons$tx_id == tx$tx_id, , drop = FALSE]; if (!nrow(e)) next
 traces[[length(traces) + 1L]] <- list(type = "scatter", mode = "lines", x = as.vector(rbind(e$start, e$end, NA)), y = as.vector(rbind(rep(i, nrow(e)), rep(i, nrow(e)), NA)), line = list(color = "#1f4e79", width = 14), connectgaps = FALSE, showlegend = FALSE, hoverinfo = "text", text = paste0("Transcript: ", tx$transcript, "<br>Exon ", e$exon_number))
 }
 list(traces = traces, names = model$transcripts$transcript)
}

build_guide_traces <- function(table) {
 if (!nrow(table)) return(list())
 colours <- c("+" = "#d95f02", "-" = "#1b9e77", "*" = "#7570b3")
 lapply(intersect(names(colours), unique(table$strand)), function(s) {
 d <- table[table$strand == s,
 d <- table[table$strand == s, , drop = FALSE]
 list(type = "scattergl", mode = "markers", name = paste(s, "strand"), x = d$cut_site, y = d$composite_score, text = paste0("Guide: ", d$guide_id, "<br>Location: ", d$seqnames, ":", d$start, "-", d$end, "<br>Strand: ", d$strand, "<br>Protospacer: ", d$protospacer, "<br>PAM: ", d$pam, "<br>GC%: ", round(d$gc_percent, 2)), hoverinfo = "text", marker = list(size = 7, opacity = 0.65, color = unname(colours[s])))
 })
}

build_table_rows <- function(table, n) {
 if (!nrow(table)) return("<tr><td colspan='11'>No guides passed the configured filters.</td></tr>")
 d <- utils::head(table[order(table$composite_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE], n)
 paste(vapply(seq_len(nrow(d)), function(i) paste0(
 "<tr><td>", escape_html(d$guide_id[i]), "</td><td>", escape_html(d$seqnames[i]), ":", escape_html(d$cut_site[i]),
 "</td><td>", escape_html(d$strand[i]), "</td><td class='sequence'>", escape_html(d$protospacer[i]),
 "</td><td class='sequence'>", escape_html(d$pam[i]), "</td><td>", escape_html(round(d$gc_percent[i], 2)),
 "</td><td>", escape_html(d$transcript[i] %||% ""), "</td><td>", escape_html(d$mature_transcript_position[i] %||% ""),
 "</td><td>", escape_html(d$transcript_region[i] %||% ""), "</td><td>", escape_html(d$exon_number[i] %||% ""),
 "</td><td>", escape_html(round(d$composite_score[i], 3)), "</td></tr>"), character(1)), collapse = "\n")
}

create_report <- function(table, model, transcript_map, pairs, config, html_file) {
 tx <- build_transcript_traces(model)
 guide_json <- jsonlite::toJSON(build_guide_traces(table), auto_unbox = TRUE, na = "null", null = "null")
 tx_json <- jsonlite::toJSON(tx$traces, auto_unbox = TRUE, na = "null", null = "null")
 map_json <- jsonlite::toJSON(transcript_map, dataframe = "rows", auto_unbox = TRUE, na = "null", null = "null")
 ticks <- jsonlite::toJSON(seq_along(tx$names), auto_unbox = TRUE)
 names_json <- jsonlite::toJSON(tx$names, auto_unbox = TRUE)
 genome <- config$genome$package %||% config$genome$source
 region <- if (!is.null(config$target$seqnames)) paste0(config$target$seqnames, ":", config$target$start, "-", config$target$end) else "gene / custom sequence"
 title <- paste(config$project_name, "sgRNA report")
 initial_transcript <- config$transcript$selected_id %||% if (length(tx$names)) tx$names[1L] else ""
 selector <- if (isTRUE(config$transcript$selector_in_report) && length(tx$names)) {
 paste0("<label for='transcriptSelect'><b>Transcript:</b></label><select id='transcriptSelect'></select>")
 } else ""
 pair_html <- if (nrow(pairs)) paste0("<h2>Paired-guide deletion candidates</h2><table><thead><tr><th>Pair</th><th>Guide 1</th><th>Guide 2</th><th>Cut 1</th><th>Cut 2</th><th>Deletion (bp)</th><th>Mean score</th></tr></thead><tbody>", paste(vapply(seq_len(min(nrow(pairs), config$output$top_guides)), function(i) paste0("<tr><td>", pairs$pair_id[i], "</td><td>", pairs$guide_1[i], "</td><td>", pairs$guide_2[i], "</td><td>", pairs$cut_1[i], "</td><td>", pairs$cut_2[i], "</td><td>", pairs$deletion_bp[i], "</td><td>", round(pairs$mean_score[i], 3), "</td></tr>"), character(1)), collapse = ""), "</tbody></table>") else ""
 html <- paste0(
 "<!doctype html><html><head><meta charset='UTF-8'><title>", escape_html(title), "</title>",
 "<script src='https://cdn.plot.ly/plotly-2.35.2.min.js'></script><style>",
 "body{font-family:Arial;margin:24px;color:#17202a}.summary{background:#f4f7f9;padding:16px;border-radius:8px;margin-bottom:20px}.plot{width:100%;height:500px;margin-bottom:28px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #d5d8dc;padding:7px;text-align:left}th{background:#1f4e79;color:white;position:sticky;top:0}tr:nth-child(even){background:#f7f9f9}.sequence{font-family:monospace}select{padding:5px;margin:8px 0 12px}</style></head><body>",
 "<h1>", escape_html(title), "</h1><div class='summary'><b>Target:</b> ", escape_html(config$target$name),
 "<br><b>Genome:</b> ", escape_html(genome), "<br><b>Region:</b> ", escape_html(region),
 "<br><b>Nuclease:</b> ", escape_html(config$nuclease$name), "<br><b>Modality:</b> ", escape_html(config$modality$type),
 "<br><b>Cell context:</b> ", escape_html(config$cell_context$name), "<br><b>Candidate guides:</b> ", nrow(table),
 "<br><b>Annotated transcripts:</b> ", if (is.null(model)) 0L else nrow(model$transcripts),
 "</div><h2>Transcript and exon structure</h2>", selector, "<div id='genePlot' class='plot'></div>",
 "<h2>Guide positions on the selected mature transcript</h2><div id='transcriptGuidePlot' class='plot'></div>",
 "<h2>Guide cut sites and scores</h2><div id='guidePlot' class='plot'></div>",
 "<h2>Top ", config$output$top_guides, " guides</h2><table><thead><tr><th>Guide</th><th>Cut site</th><th>Strand</th><th>Protospacer</th><th>PAM</th><th>GC %</th><th>Transcript</th><th>Mature transcript position</th><th>Region</th><th>Exon</th><th>Composite score</th></tr></thead><tbody>",
 build_table_rows(table, config$output$top_guides), "</tbody></table>", pair_html,
 "<script>const geneData=", tx_json, ";const guideData=", guide_json, ";const transcriptMap=", map_json, ";const transcriptNames=", names_json,
 ";const plotConfig={responsive:true,displaylogo:false,scrollZoom:true};",
 "const geneLayout={title:'Transcripts and exons',xaxis:{title:'Genomic position',rangeslider:{visible:true}},yaxis:{title:'Transcript',tickmode:'array',tickvals:", ticks, ",ticktext:", names_json, "},margin:{l:180,r:30,t:60,b:70}};Plotly.newPlot('genePlot',geneData,geneLayout,plotConfig);",
 "Plotly.newPlot('guidePlot',guideData,{title:'Candidate guides',xaxis:{title:'Cut-site position'},yaxis:{title:'Composite score'},margin:{l:80,r:30,t:60,b:70}},plotConfig);",
 "function drawTranscriptGuides(tx){const rows=transcriptMap.filter(x=>x.transcript===tx&&x.transcript_position!==null);const strands=['+','-','*'];const colors={'+':'#d95f02','-':'#1b9e77','*':'#7570b3'};const traces=strands.map(s=>{const r=rows.filter(x=>x.strand===s);return {type:'scattergl',mode:'markers',name:s+' strand',x:r.map(x=>x.transcript_position),y:r.map(x=>x.composite_score),text:r.map(x=>'Guide: '+x.guide_id+'<br>Exon: '+x.exon_number+'<br>Region: '+x.transcript_region+'<br>Position: '+x.transcript_position+'<br>Distance to TSS: '+x.distance_to_tss),hoverinfo:'text',marker:{size:8,color:colors[s]}}});Plotly.react('transcriptGuidePlot',traces,{title:'Cleavage positions on '+tx,xaxis:{title:'Mature transcript position (1-based)'},yaxis:{title:'Composite score'},margin:{l:80,r:30,t:60,b:70}},plotConfig);const txIndex=transcriptNames.indexOf(tx)+1;const genomicTraces=strands.map(s=>{const r=transcriptMap.filter(x=>x.transcript===tx&&x.strand===s&&x.cut_site!==null);return {type:'scattergl',mode:'markers',name:s+' sgRNA cuts',x:r.map(x=>x.cut_site),y:r.map(x=>txIndex),text:r.map(x=>'Guide: '+x.guide_id+'<br>Genomic cut: '+x.seqnames+':'+x.cut_site+'<br>Mature transcript position: '+x.transcript_position+'<br>Score: '+x.composite_score),hoverinfo:'text',marker:{size:9,color:colors[s],line:{color:'#17202a',width:0.5}}};});Plotly.react('genePlot',geneData.concat(genomicTraces),geneLayout,plotConfig);} ",
 "const transcriptSelect=document.getElementById('transcriptSelect');if(transcriptSelect){const txs=[...new Set(transcriptMap.map(x=>x.transcript))].filter(Boolean);txs.forEach(tx=>{const o=document.createElement('option');o.value=tx;o.textContent=tx;transcriptSelect.appendChild(o);});const requested=", jsonlite::toJSON(initial_transcript, auto_unbox=TRUE), ";const initial=txs.includes(requested)?requested:(transcriptSelect.options.length?transcriptSelect.options[0].value:null);if(initial){transcriptSelect.value=initial;drawTranscriptGuides(initial);}transcriptSelect.addEventListener('change',e=>drawTranscriptGuides(e.target.value));}else if(transcriptMap.length){drawTranscriptGuides(transcriptMap[0].transcript);} </script></body></html>"
 )
 writeLines(html, html_file, useBytes = TRUE); invisible(html_file)
}

run_pipeline <- function(config = default_config()) {
 config <- validate_config(config); install_required_packages(config)
 bsgenome <- resolve_bsgenome(config); nuclease <- resolve_nuclease(config)
 guides <- design_guides(config, nuclease, bsgenome)
 if (!is.null(config$cell_context$vcf) && config$mode == "design") {
 guides <- try_add(guides, crisprDesign::addSNPAnnotation,
 vcf = config$cell_context$vcf, label = "Cell-context VCF annotation")
 }
 if (isTRUE(config$off_target$enabled) && config$mode == "design") {
 if (is.null(config$off_target$bowtie_index)) {
 warning("off_target$enabled is TRUE but no bowtie_index was supplied; off-target search skipped.", call. = FALSE)
 } else {
 guides <- try_add(guides, crisprDesign::addSpacerAlignments,
 aligner_index = config$off_target$bowtie_index,
 n_mismatches = config$off_target$mismatches,
 all_alignments = config$off_target$all_alignments,
 label = "Off-target alignment")
 guides <- try_add(guides, crisprDesign::addOffTargetScores,
 label = "Off-target scoring")
 }
 }
 if (config$mode == "design" && isTRUE(config$scoring$add_sequence_features)) guides <- try_add(guides, crisprDesign::addSequenceFeatures, label = "Sequence-feature scoring")
 if (config$mode == "design" && isTRUE(config$scoring$add_on_target)) {
 guides <- try_add(guides, crisprDesign::addOnTargetScores,
 methods = config$scoring$methods %||% c("ruleset1", "crisprater"),
 label = "On-target scoring")
 }
 model <- if (isTRUE(config$annotation$enabled)) get_gene_model(config) else NULL
 mapped <- annotate_transcript_positions(normalise_guides(guides, config), model, config)
 table <- annotate_modality(mapped$table, config)
 table <- annotate_cloning(table, config)
 table <- filter_guides(table, config)
 pairs <- build_guide_pairs(table, config)
 dir.create(config$output$directory, recursive = TRUE, showWarnings = FALSE)
 prefix <- file.path(config$output$directory, config$output$prefix)
 paths <- list(
 rds = paste0(prefix, ".rds"), csv = paste0(prefix, ".csv"), html = paste0(prefix, ".html"),
 bed = paste0(prefix, ".bed"), pairs = paste0(prefix, "_pairs.csv"),
 transcript_map = paste0(prefix, "_transcript_map.csv")
 )
 if (isTRUE(config$output$write_rds)) saveRDS(table, paths$rds)
 if (isTRUE(config$output$write_csv)) utils::write.csv(table, paths$csv, row.names = FALSE)
 if (isTRUE(config$output$write_bed)) write_bed(table, paths$bed)
 if (isTRUE(config$output$write_pairs)) utils::write.csv(pairs, paths$pairs, row.names = FALSE)
 if (isTRUE(config$output$write_transcript_map)) utils::write.csv(mapped$map, paths$transcript_map, row.names = FALSE)
 if (isTRUE(config$output$write_html)) create_report(table, model, mapped$map, pairs, config, paths$html)
 message("Designed ", nrow(table), " guides."); if (isTRUE(config$output$write_html)) message("Interactive report: ", paths$html)
 invisible(list(config = config, guides = table, guide_set = guides, gene_model = model, transcript_map = mapped$map, pairs = pairs, paths = paths))
}

if (identical(environment(), globalenv())) {
 args <- commandArgs(trailingOnly = TRUE)
 cfg <- default_config()
 if (length(args)) cfg <- load_config_file(args[1L], cfg)
 # source("...") in an interactive R session now produces the report
 # instead of only defining functions with no visible result.
 run_pipeline(cfg)
}
