Trumpet_report <- function(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
                           condition1, condition2, GENE_ANNO_GTF = NA, TXDB = NA, sample_size = NA, 
                           GENOME = NA, UCSC_TABLE_NAME = "knownGene", 
                           OUTPUT_DIR = NA) {
    outparam_dir <- tempdir()
    if (is.na(OUTPUT_DIR)) {
      OUTPUT_DIR <- getwd()
      trumpet <- system.file("extdata", "Trumpet_report.Rmd", package = "Trumpet")
      if (suppressWarnings(is.na(GENOME) & is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
        
        stop("Please give the annotation gene in GTF format or TXDB file or download the knownGene from the UCSC")
      }
      if (suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME)) & 
                           is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
        if (is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENOME, UCSC_TABLE_NAME, GENE_ANNO_GTF, 
               TXDB, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                           "parameter.Rdata", sep = "/"))
        }
        if (!is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENOME, UCSC_TABLE_NAME, GENE_ANNO_GTF, 
               TXDB, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                           "parameter.Rdata", sep = "/"))
        }
      }
      if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB))) {
        if (is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        if (!is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
      }
      if (suppressWarnings(!is.na(TXDB))) {
        if (is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, TXDB, GENE_ANNO_GTF, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        if (!is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, TXDB, GENE_ANNO_GTF, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        
      }
      render(trumpet, output_format = "html_document", output_dir = OUTPUT_DIR)
    }
    if (!is.na(OUTPUT_DIR)) {
      trumpet <- system.file("extdata", "Trumpet_report.Rmd", package = "Trumpet")
      if (suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME)) & 
                           is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
        if (is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        if (!is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
      }
      if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB))) {
        if (is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        if (!is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
      }
      if (suppressWarnings(!is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
        if (is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        if (!is.na(sample_size)) {
          save(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
               condition1, condition2, GENE_ANNO_GTF, TXDB, GENOME, 
               UCSC_TABLE_NAME, sample_size, OUTPUT_DIR, file = paste(outparam_dir, 
                                                                      "parameter.Rdata", sep = "/"))
        }
        
      }
      render(trumpet, output_format = "html_document", output_dir = OUTPUT_DIR)
    }
}
