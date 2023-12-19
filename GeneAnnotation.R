gene_annotation <- function(Model_name,gff_FileName, snp_FileName, Gwas.output_FileName,range){
    gff = read.table(gff_FileName,header = TRUE , sep="\t")
    snp = read.csv(snp_FileName)
    result <- list()
    dirname = paste0(Model_name,"_genes_output")
    dir.create(dirname)
    
    for(i in 1:nrow(snp)){
        snp_num = snp$Position[i]
        snp_name = snp$SNP[i]
        snp_chr = snp$Chromosome[i]
        vector_range = gff$start[(gff$start > snp_num - range) & (gff$start < snp_num + range)]
        ranged_genes = subset(gff,start %in% vector_range)
        ranged_genes = ranged_genes[ranged_genes$chr == snp_chr,]
        
        df_name <- paste0("df_", snp_name, "_",Model_name)
        result[[df_name]] <- ranged_genes
        
        write.csv(ranged_genes,paste0(dirname,"/",paste0(df_name,".csv")))
        
    }
    
    #manheten plot Gwas
    
    df = read.csv(Gwas.output_FileName)
    library(qqman)
    
    df$BP <- as.numeric(rownames(df))
    df$P <- df$P.value
    df$CHR <- as.numeric(gsub("chr(\\d+)H", "\\1", df$Chr))
    
    pdf(paste0(dirname,"/",Model_name , "_Manhattan.pdf"))
    manhattan(df, annotatePval = 0.01)
    title(main = Model_name)
    dev.off()
    
    return(result)

}
