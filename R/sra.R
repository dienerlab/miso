# Helpers for creating SRA submissions

instruments <- list(
    ILLUMINA = c("HiSeq X Five",
                 "HiSeq X Ten",
                 "Illumina Genome Analyzer",
                 "Illumina Genome Analyzer II",
                 "Illumina Genome Analyzer IIx",
                 "Illumina HiScanSQ",
                 "Illumina HiSeq 1000",
                 "Illumina HiSeq 1500",
                 "Illumina HiSeq 2000",
                 "Illumina HiSeq 2500",
                 "Illumina HiSeq 3000",
                 "Illumina HiSeq 4000",
                 "Illumina iSeq 100",
                 "Illumina NovaSeq 6000",
                 "Illumina MiniSeq",
                 "Illumina MiSeq",
                 "NextSeq 500",
                 "NextSeq 550"),
    HELICOS = c("Helicos HeliScope"),
    ABI_SOLID = c("AB 5500 Genetic Analyzer",
                  "AB 5500xl Genetic Analyzer",
                  "AB 5500x-Wl Genetic Analyzer",
                  "AB SOLiD 3 Plus System",
                  "AB SOLiD 4 System",
                  "AB SOLiD 4hq System",
                  "AB SOLiD PI System",
                  "AB SOLiD System",
                  "AB SOLiD System 2.0",
                  "AB SOLiD System 3.0"),
    COMPLETE_GENOMICS = c("Complete Genomics"),
    PACBIO_SMRT = c("PacBio RS",
                    "PacBio RS II",
                    "PacBio Sequel"),
    ION_TORRENT = c("Ion Torrent PGM",
                    "Ion Torrent Proton",
                    "Ion Torrent S5 XL",
                    "Ion Torrent S5"),
    CAPILLARY = c("AB 310 Genetic Analyzer",
                  "AB 3130 Genetic Analyzer",
                  "AB 3130xL Genetic Analyzer",
                  "AB 3500 Genetic Analyzer",
                  "AB 3500xL Genetic Analyzer",
                  "AB 3730 Genetic Analyzer",
                  "AB 3730xL Genetic Analyzer"),
    OXFORD_NANOPORE = c("GridION", "MinION", "PromethION"),
    BGISEQ = c("BGISEQ-500")
)

gut_envo <- c(
    env_broad_scale = paste0("environment associated with an animal part ",
                             "or small animal [ENVO:01001055]"),
    env_local_scale = "intestine environment [ENVO:2100002]",
    env_medium = "fecal material [ENVO:00002003]"
)

skin_envo <- c(
    env_broad_scale = paste0("environment associated with an animal part ",
                             "or small animal [ENVO:01001055]"),
    env_local_scale = "integumental system environment [ENVO:2100004]",
    env_medium = "skin environment [ENVO:2100003]"
)

water_biofilm_envo <- c(
    env_broad_scale = "freshwater environment [ENVO:01000306]",
    env_local_scale = paste0("environment determined by a biofilm ",
                             "on a non-saline surface [ENVO:01001051]"),
    env_medium = "biofilm material [ENVO:01000156]"
)

in_vitro_envo <- c(
    env_broad_scale = "organic material [ENVO:01000155]",
    env_local_scale = "mock community culture [ENVO:01001059]",
    env_medium = "cell culturing unit [ENVO:01001819]"
)

usage_16S <- paste0(
    "You are now ready for submission. ",
    "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
    "log in and click on `New submission`. ",
    "Fill in the general data for your project in steps ",
    "1 through 3. In step 4 ",
    "choose `Genome, metagenome or marker sequences ",
    "(MIxS compliant)` and `Survey-related Marker ",
    "Sequences MIMARKS`. In step 5 and 6 ",
    "upload the respective files in %s. In step 7 you ",
    "can directly upload the `*.tar.gz` submission ",
    "package. Just click on `continue` another time to ",
    "have the archive unpacked as indicated."
)

usage_mx <- paste0(
    "You are now ready for submission. ",
    "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
    "log in and click on `New submission`. ",
    "Fill in the general data for your project in steps ",
    "1 through 3. In step 4 ",
    "choose `Environmental/Metagenome ",
    "Genomic Sequences MIMS` with the appropriate environment. In step 5 and 6 ",
    "upload the respective files in %s. In step 7 you ",
    "can directly upload the `*.tar.gz` submission ",
    "package. Just click on `continue` another time to ",
    "have the archive unpacked as indicated."
)

#' The available presets for SRA submissions. They all map to MIMARKS or MIMS
#' packages
#'
#' @export
sra_presets <- list(
    `human gut 16S` = c(
        gut_envo,
        organism = "human gut metagenome",
        host = "Homo sapiens",
        library_strategy = "AMPLICON",
        library_source = "GENOMIC",
        library_selection = "PCR",
        filetype = "fastq",
        usage = usage_16S
    ),
    `mouse gut 16S` = c(
        gut_envo,
        organism = "mouse gut metagenome",
        host = "Mus musculus",
        library_strategy = "AMPLICON",
        library_source = "GENOMIC",
        library_selection = "PCR",
        filetype = "fastq",
        usage = usage_16S
    ),
    `human gut metagenome` = c(
        gut_envo,
        organism = "human gut metagenome",
        host = "Homo sapiens",
        library_strategy = "WGS",
        library_source = "METAGENOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = usage_mx
    ),
    `human skin 16S` = c(
        skin_envo,
        organism = "human skin metagenome",
        host = "Homo sapiens",
        library_strategy = "AMPLICON",
        library_source = "GENOMIC",
        library_selection = "PCR",
        filetype = "fastq",
        usage = usage_16S
    ),
    `human skin metagenome` = c(
        skin_envo,
        organism = "human skin metagenome",
        host = "Homo sapiens",
        library_strategy = "WGS",
        library_source = "METAGENOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = usage_mx
    ),
    `mouse gut metagenome` = c(
        gut_envo,
        organism = "mouse gut metagenome",
        host = "Mus musculus",
        library_strategy = "WGS",
        library_source = "METAGENOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = usage_mx
    ),
    `human gut RNA-Seq` = c(
        gut_envo,
        organism = "human gut metagenome",
        host = "Homo sapiens",
        library_strategy = "RNA-Seq",
        library_source = "METATRANSCRIPTOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = usage_mx
    ),
    `mouse gut RNA-Seq` = c(
        gut_envo,
        organism = "mouse gut metagenome",
        host = "Mus musculus",
        library_strategy = "RNA-Seq",
        library_source = "METATRANSCRIPTOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = usage_mx
    ),
     `aquatic biofilm 16S` = c(
        water_biofilm_envo,
        organism = "environmental metagenome",
        depth = "0 m",
        elevation = "840 m",  # global average as placeholder
        library_strategy = "AMPLICON",
        library_source = "GENOMIC",
        library_selection = "PCR",
        filetype = "fastq",
        usage = usage_16S
    ),
    `in vitro metagenome` = c(
        in_vitro_envo,
        organism = "mixed culture metagenome",
        library_strategy = "WGS",
        library_source = "METAGENOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = usage_mx
    )
)

validate <- function(config) {
    if (is.null(config$metadata)) {
        stop("Need to specify a metadata files in config :(")
    }
    if (is.null(config$title)) {
        stop(paste0("Please specify a title in the form ",
                    "`{methodology} of {organism}: {sample info}`"))
    }
    if (!config$preset %in% names(sra_presets)) {
        stop(sprintf("Not a supported presets. Allowed are: %s",
                     paste0(names(sra_presets), collapse = ", ")))
    }
    if (!toupper(config$platform) %in% names(instruments)) {
        stop(sprintf("Not a recognized platform, please pick one of: %s",
             paste0(names(instruments), collapse = ",")))
    }
    if (!config$instrument_model %in%
        instruments[[toupper(config$platform)]]) {
        stop(sprintf("Not a recognized model, please pick one of: %s",
             paste0(instruments[[toupper(config$platform)]], collapse = ",")))
    }
}

#' Build a configuration for the SRA submission workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the transcript counting
#'  workflow.
#' @export
#' @examples
#'  config <- config_sra(date_col = "collection_date")
config_sra <- config_builder(list(
    metadata = NULL,
    id_col = "id",
    date_col = "date",
    country = "USA",
    preset = "human gut 16S",
    out_dir = "sra",
    title = NULL,
    platform = "ILLUMINA",
    instrument_model = "Illumina MiSeq",
    bioproject = NULL,
    latitude = 42.36,
    longitude = -71.0941,
    make_package = TRUE
))


#' Prepare submission files for the NCBI sequence read archive (SRA).
#'
#' @param object A files data table as returned by \code{\link{find_read_files}}.
#' @param ... Configuration arguments or a configuration as generated by
#'  \code{\link{config_sra}}.
#' @return A list containing the used alignments and the transcript counts in
#'  `counts`.
#'
#' @export
#' @importFrom utils tar
#' @importFrom stringr str_replace_all
sra_submission <- function(object, ...) {
    files <- get_files(object)
    config <- config_parser(list(...), config_sra)
    validate(config)
    if (!dir.exists(config$out_dir)) {
        dir.create(config$out_dir, recursive = TRUE)
    }
    meta <- as.data.table(config$metadata)
    preset <- sra_presets[[config$preset]]
    files <- copy(files)
    setkey(files, id)
    if (!all(meta[[config$id_col]] %in% files$id)) {
        stop("Some ids in the metadata have no matching files :(")
    }
    if (!all(files$id %in% meta[[config$id_col]])) {
        flog.info(paste0("Some file ids are not mentioned in the metadata. ",
                         "Will assume that is intended and omit those."))
    }
    files <- files[config$metadata[[config$id_col]]]
    date <- as.Date(meta[[config$date_col]])
    date <- format(date, "%Y-%m-%d")
    sample_data <- data.table(
        sample_name = files$id,
        organism = preset["organism"],
        collection_date = date,
        env_broad_scale = preset["env_broad_scale"],
        env_local_scale = preset["env_local_scale"],
        env_medium = preset["env_medium"],
        host = preset["host"],
        geo_loc_name = config$country,
        lat_lon = paste(round(abs(config$latitude), 4),
                        ifelse(config$latitude >= 0, "N", "S"),
                        round(abs(config$longitude), 4),
                        ifelse(config$longitude >= 0, "E", "W")),
        source_material_id = files$id
    )

    if ("host" %in% names(preset)) {
        sample_data$host <- preset["host"]
    } else if ("depth" %in% names(preset)) {
        sample_data$depth <- preset["depth"]
        sample_data$elevation <- preset["elevation"]
    }
    additional <- names(meta)[!names(meta) %in%
                              c(config$id_col, config$date_col)]
    additional <- meta[, additional, with = FALSE]
    names(additional) <- tolower(str_replace_all(names(additional),
                                 "[^A-Za-z0-9_]", "_"))
    sample_data <- cbind(sample_data, additional)

    sra_metadata <- data.table(
        sample_name = files$id,
        library_ID = files$id,
        filename = basename(files$forward),
        title = config$title,
        library_strategy = preset["library_strategy"],
        library_source = preset["library_source"],
        library_selection = preset["library_selection"],
        filetype = preset["filetype"],
        library_layout = ifelse("reverse" %in% names(files),
                                "paired", "single"),
        platform = toupper(config$platform),
        instrument_model = config$instrument_model,
        design_description = "see manuscript for methods"
    )
    if (!is.null(config$bioproject)) {
        sample_data[, "bioproject_accession" := config$bioproject]
        sra_metadata[, "bioproject_accession" := config$bioproject]
    }
    if ("reverse" %in% names(files)) {
        sra_metadata[, "filename2" := basename(files$reverse)]
    }

    packages <- list()
    upload <- list()
    if (nrow(sample_data) < 1000) {
        upload[[1]] <- files$forward
        if ("reverse" %in% names(files)) {
            upload[[1]] <- c(upload[[1]], files$reverse)
        }
        if (config$make_package) {
            if (Sys.getenv("tar") == "") {
                tar <- Sys.getenv("TAR")
            } else {
                tar <- Sys.getenv("tar")
            }
            packages[[1]] <- file.path(config$out_dir, "sra_files.tar.gz")
            flog.info("Packing submission files to %s.",
                    file.path(config$out_dir, "sra_files.tar.gz"))
            ret <- tar(file.path(config$out_dir, "sra_files.tar.gz"),
                files = upload[[1]],
                compression = "gzip", tar = tar)
            if (ret != 0) {
                stop("tar command failed")
            }
        }
        flog.info("Writing biosample attributes to %s.",
                file.path(config$out_dir, "05_biosample_attributes.tsv"))
        fwrite(sample_data, sep = "\t", file.path(config$out_dir,
                                                "05_biosample_attributes.tsv"))
        flog.info("Writing file metadata to %s.",
                file.path(config$out_dir, "06_sra_metadata.tsv"))
        fwrite(sra_metadata, sep = "\t",
            file.path(config$out_dir, "06_sra_metadata.tsv"))
    } else {
        start <- 1
        # distribute samples evenly over submission
        nsub <- ceiling(nrow(sample_data) / 999)
        step <- ceiling(nrow(sample_data) / nsub)
        flog.info(paste0(
            "SRA does not allow submissions with more than 1,000 samples. ",
            "I will split your submission into ", nsub, " chunks, which ",
            "you will have to submit separately. Using the instructions below."
        ))
        chunk <- 1
        while (start < nrow(sample_data)) {
            idx <- start : min((start + step - 1), nrow(sample_data))
            upload[[chunk]] <- files$forward[idx]
            if ("reverse" %in% names(files)) {
                upload[[chunk]] <- c(upload[[chunk]], files$reverse[idx])
            }
            if (config$make_package) {
                if (Sys.getenv("tar") == "") {
                    tar <- Sys.getenv("TAR")
                } else {
                    tar <- Sys.getenv("tar")
                }
                out <- sprintf("sra_files_%d.tar.gz", chunk)
                packages[[chunk]] <- file.path(config$out_dir, out)
                flog.info("Packing submission files for chunk %d to %s.",
                        chunk, file.path(config$out_dir, out))
                ret <- tar(file.path(config$out_dir, out),
                    files = upload[[chunk]],
                    compression = "gzip", tar = tar)
                if (ret != 0) {
                    stop("tar command failed")
                }
            }
            satts <- sprintf("05_biosample_attributes_%d.tsv", chunk)
            flog.info("Writing biosample attributes for chunk %d to %s.",
                chunk, file.path(config$out_dir, satts))
            fwrite(sample_data[idx], sep = "\t",
                file.path(config$out_dir, satts))
            smeta <- sprintf("06_sra_metadata_%d.tsv", chunk)
            flog.info("Writing file metadata for chunk %d to %s.",
                    chunk, file.path(config$out_dir, smeta))
            fwrite(sra_metadata[idx], sep = "\t",
                file.path(config$out_dir, smeta))

            chunk <- chunk + 1
            start <- start + step
        }
    }
    flog.info(sprintf(preset["usage"], config$out_dir))
    artifact <- list(
        files = files,
        biosample_attributes = sample_data,
        sra_metadata = sra_metadata
    )
    if (config$make_package) {
        artifact$upload <- packages
    } else {
        artifact$upload <- upload
    }
    return(artifact)
}
