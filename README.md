# Gene_annotation
this library tryies to solve the problem of gene annotation or after GWAS analisys wich there areent any good library who does that.

## what the function do
the function outputs a folder that contain tables significant genes in a window of nucleotides that recived from the user and manhetten plot with the SNP annotation
### how the folder looks like:
![image](https://github.com/talamar49/Gene_annotation/assets/114323965/5dd1dcd8-be51-4288-b2ef-2514c4800ada)
![image](https://github.com/talamar49/Gene_annotation/assets/114323965/07f6820b-8576-4973-b543-98845677a123)

### how the output table looks like:
the snps tables:

![image](https://github.com/talamar49/Gene_annotation/assets/114323965/c8784e52-4d3e-4b3e-a1d7-0e286f93c5b1)

the Manhattan plot pdf:

![image](https://github.com/talamar49/Gene_annotation/assets/114323965/80647f7b-a710-4fa8-89e6-c2e46e01a851)


## what argument the function needs

! keep in mind to run this script inside your Working directory with all the files mentioned in it

the function needs:

1. model name: how you want your folder/snp subname/plot title, to be called.
2. gff_FileName: a string of the gff (gene annotation file)
3. snp_FileName: a string of significant snps (come with the output from the GWAS files)
4. Gwas.output_Filename: the name of the Gwas output (for manhattan plot)
5. range: the window of nucleotide to search for genes around the given snp

this is how the file looks in the folder.

![image](https://github.com/talamar49/Gene_annotation/assets/114323965/93d9ec20-5ca8-43e7-9e73-b83a5dbf7359)

the first on is the gff converted to txt file this is how it looks like the output csv file look at "how the output table looks like"

the second file is the GWAS result , the results came from gapit (https://www.maizegenetics.net/gapit) with the model BLINK and it looks like this:

![image](https://github.com/talamar49/Gene_annotation/assets/114323965/201a4053-5940-41f0-8e4e-9b304e85399c)

the third file is the significant snp csv and it looks like this:

![image](https://github.com/talamar49/Gene_annotation/assets/114323965/7dbda646-a94a-4816-942a-7b216bcc54bb)

## the function itself
```R
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
```

## function use:

general:
```R
gene_annotation(Model_name ,gff_FileName , snp_FileName , Gwas.output_FileName ,range)
```
example:
```R
gene_annotation("BLINK_beta", gff_FileName , snp_beta , Gwas_beta , 1000000)
```

## what the algorithem:
basicly the function takes the snp position from the significant snp file and searches for genes in the txt(gff converted) file in the range specified (in the example above 1M window) and saves the dataframe into a folder.
