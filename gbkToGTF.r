# Converting a genebank file to a STAR&RSEM ready gtf file
# Benjamin Saintpierre
# Plateforme GENOM'IC


# reading a genebank
# BiocManager::install("genbankr")

suppressPackageStartupMessages(library(genbankr))

# Lecture d'un fichier genebank issu du NCBI
# Reading a NCBI genebank file
gb <- readGenBank("/path/to/gbk/file.gbk")

# Fonction pour lire et extraire les informations d'un dataframe issu de genbankr::readGenBank()
# Function to read and extract informations from a genbankr::readGenBank() dataframe
convert_genesGb_gtf <- function(ligne_genebank){
  
  # Vérification de la présence d'un gene_id. Dans le cas contraire, récuperer le label.
  # Attention, toujours vérifier en amont qu'au moins un de ces deux éléments est disponible pour le gène
  # Dans le cas contraire, identifier une autre colonne non-vide, susceptible de se subsituer au gene_id.
  
  # Getting the gene_id value, or the label if the gene_id is not available
  # Please be sure that one of these values is not empty
  # Otherwise you will need to modify the label column to a new column get a right value 
  
  if(is.na(ligne_genebank$gene_id)){
    geneID <- ligne_genebank$label
  }else{
    geneID <- ligne_genebank$gene_id
  }
  
  # Récupération des éléments necessaires à la génération d'un gtf
  # Basé sur la composition obligatoire d'une ligne d'un fichier gtf.
  
  # Getting all gtf mandatory elements
  geneName <- ligne_genebank$gene
  gene_biotype <- ligne_genebank$type
  gene_start <- ligne_genebank$start
  gene_end <- ligne_genebank$end
  gene_strand <- ligne_genebank$strand
  
  # La description est importante, doit etre composé d'éléments séparés par des ;
  # chaque éléement doit avoir la construction suivante: nom_de_l_element "valeur"
  # l'espace et les guillemets sont indispensables
  # Notre pipeline utilise STAR et RSEM, donc éléments obligatoires: gene_id & transcript_id
  
  # Description column matter a lot, a muste be composed as:
  # element_name "value";
  # Important: space and quotation marks must be there
  # gene_id & transcript_id have to be part of the informations for our STAR + RSEM pipeline
  description <- paste("gene_id \"",geneID,"\"; transcript_id \"",geneID,"\"; gene_name \"",geneName,"\"; gene_source \"","custom ","\"; gene_biotype \"",gene_biotype,"\"",sep="")
  
  # Retourne un vecteur pour chaque vlaure.
  # Attention!! Penser à modifier la première colonne pour la faire correspondre au nom du chromosome
  # Si plusieurs valeurs sont possibles, changer cette valeur pour y intégrer une variable, sur l'exemple des éléments ci-dessus
  
  # Return a vector with each information
  # Be careful, here the chromosom name is written, not a variable
  # If many values are available, please make a new line and variable, based on the how we get element here
  return(c("NC_007793.1","custom","exon",gene_start,gene_end,".",gene_strand,".",description))
}

# Conversion en dataframe nécessaire pour l'utilisation de la fonction apply
# From GenBankRecord object to dataframe
genes_gb <- as.data.frame(genes(gb))
# La transposition est nécessaire ici
# Transposition is needed
new_gtf_file <- t(apply(X = genes_gb,1,convert_genesGb_gtf))

# Sauvegarde d'un fichier gtf, utilisable avec STAR et RSEM
# Saving the STAR & RSEM ready gtf file
write.table(x = new_gtf_file,file = "genes.gtf",sep="\t", quote = FALSE, col.names = FALSE,row.names = FALSE)




# Même fonction mais pour les annotations non "gene"
# même remarques, et même fonctionnement que précédemment

# Same function for not-gene annotations

convert_otherFeaturesGb_gtf <- function(ligne_genebank){
  
  if(is.na(ligne_genebank$locus_tag)){
    geneID <- ligne_genebank$note
  }else{
    geneID <- ligne_genebank$locus_tag
  }
  geneName <- ligne_genebank$gene
  gene_biotype <- ligne_genebank$type
  gene_start <- ligne_genebank$start
  gene_end <- ligne_genebank$end
  gene_strand <- ligne_genebank$strand
  
  description <- paste("gene_id ",geneID,"; gene_name ",geneName,"; gene_source ","custom ","; gene_biotype ",gene_biotype,sep="")
  
  return(c("NC_007793.1","custom","exon",gene_start,gene_end,".",gene_strand,".",description))
}

otherFeatures_gb <- as.data.frame(otherFeatures(gb))
new_gtf_file <- t(apply(X = otherFeatures_gb,1,convert_otherFeaturesGb_gtf))
write.table(x = new_gtf_file,file = "otherFeatures_gb.gtf",sep="\t", quote = FALSE, col.names = FALSE,row.names = FALSE)
