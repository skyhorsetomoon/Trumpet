Trumpet_report<-function(IP_BAM,
                         Input_BAM,
                         contrast_IP_BAM,
                         contrast_Input_BAM,
                         condition1,
                         condition2,
                         GENE_ANNO_GTF=NA,
                         TXDB = NA,
                         sample_size=NA,
                         GENOME = NA,
                         UCSC_TABLE_NAME = "knownGene",
                         get_read=FALSE,
                         INPUT_DIR=NA,
                         OUTPUT_DIR=NA){
  get_read<- as.numeric(get_read)
  if(get_read==0&is.na(INPUT_DIR))
  {
    outparam_dir<-tempdir()
    if (is.na(OUTPUT_DIR)) {
      OUTPUT_DIR<-getwd()
      trumpet<-system.file("extdata", "Trumpet_report.Rmd", package="Trumpet")
      if ( suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME))& is.na(TXDB)&is.na(GENE_ANNO_GTF))) {
        if(is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENOME,UCSC_TABLE_NAME,GENE_ANNO_GTF,TXDB,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        if(!is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENOME,UCSC_TABLE_NAME,GENE_ANNO_GTF,TXDB,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
      }
      if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB)) ){
        if(is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        if(!is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
      }
      if ( suppressWarnings(!is.na(TXDB)) )
      {
        if(is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,TXDB,GENE_ANNO_GTF,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        if(!is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,TXDB,GENE_ANNO_GTF,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        
      }
      render( trumpet,output_format ="html_document",output_dir=OUTPUT_DIR)
    }
    if (!is.na(OUTPUT_DIR)){
      trumpet<-system.file("extdata", "Trumpet_report.Rmd", package="Trumpet")
      if ( suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME))& is.na(TXDB)&is.na(GENE_ANNO_GTF))) {
        if(is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        if(!is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
      }
      if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB)) ){
        if(is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        if(!is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
      }
      if ( suppressWarnings(!is.na(TXDB) &is.na(GENE_ANNO_GTF) ) )
      {
        if(is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        if(!is.na(sample_size))
        {
          save(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM,condition1,condition2,GENE_ANNO_GTF,TXDB,GENOME,UCSC_TABLE_NAME,sample_size,OUTPUT_DIR,file=paste(outparam_dir,'parameter.Rdata',sep='/'))
        }
        
      }
      render( trumpet,output_format ="html_document",output_dir=OUTPUT_DIR)
    }
  }  
  if(get_read==1&(!is.na(INPUT_DIR)))
  {
    if (is.na(OUTPUT_DIR)) {
      trumpet<-system.file("extdata", "Trumpet_report2.Rmd", package="Trumpet")
      OUTPUT_DIR<-setwd(INPUT_DIR)
      render( trumpet,output_format ="html_document",output_dir=OUTPUT_DIR)
    }
    if (!is.na(OUTPUT_DIR)){
      trumpet<-system.file("extdata", "Trumpet_report2.Rmd", package="Trumpet")
      render( trumpet,output_format ="html_document",output_dir=OUTPUT_DIR)
    }
  }
}
