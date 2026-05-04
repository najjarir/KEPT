remove_problem_genes <- function(genes) {
  sapply(strsplit(genes, ';'), function(items) {
    filtered <- items[!grepl("[%?$]", items)]
    cleaned <- gsub("[*^]", "", filtered)
    cleaned <- cleaned[nzchar(cleaned)]
    out <- paste(cleaned, collapse = ';')

    if (out == '' || out == '-') {
      return(NA_character_)
    }
    return(out)
  })
}

exclude_genes <- function(genes, to_exclude) {
  sapply(strsplit(genes, ';'), function(items) {
    filtered <- items[!items %in% to_exclude]
    result <- paste(filtered, collapse = ';')

    if (result == '' || result == 'NA') {
      return(NA_character_)
    }

    return(result)
  })
}

create_gene_cluster <- function(genes, cluster_genes, cluster_name) {
  sapply(strsplit(genes, ';'), function(items) {
    all_present <- all(cluster_genes %in% items)
    filtered <- items[!items %in% cluster_genes]

    if (all_present) {
      filtered <- c(filtered, cluster_name)
    }

    result <- paste(filtered, collapse = ';')

    if (result == '' || result == 'NA') {
      return(NA_character_)
    }

    return(result)
  })
}


mark_truncated_genes <- function(genes) {
  genes <- sapply(genes, function(gene) {
      if (grepl("-[0-9]+%", gene)) {
        return("truncated")
      } else {
        return(gene)
      }

    result <- paste(filtered, collapse = ';')

    if (result == '' || result == 'NA') {
      return(NA_character_)
    }

    return(result)
  })
}



## Function to extract unique genes from a column
#extract_unique_genes <- function(column) {
  ## Split all entries by semicolon, unlist, get unique, sort
  #unique_genes <- column |>
    #na.omit() |>
    #strsplit(";") |>
    #unlist() |>
    #unique() |>
    #sort()
  
  #return(unique_genes)
#}


