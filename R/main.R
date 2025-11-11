
#' Encode DNA sequences into nucleotide position features
#'
#' This function converts a vector of DNA sequences (strings of "A", "T", "C", "G")
#' into a data frame representation suitable for machine learning input.
#'
#' @param dna_strings A character vector, where each element is a DNA sequence of equal length.
#'
#' @return A data frame where each column corresponds to a nucleotide position (e.g., nt_pos1, nt_pos2, ...),
#' and each cell contains a factor with levels "A", "T", "C", and "G".
#'
#' @examples
#' dna_encoding(c("ATCGT", "TGCAT"))
#'
#' @import randomForest
#' @export
#' @importFrom stats predict
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A sites for multiple RNA samples
#'
#' This function uses a trained machine learning model to predict m6A modification
#' probabilities and statuses for multiple RNA sequence entries provided in a data frame.
#'
#' @param ml_fit A trained machine learning model (e.g., from \pkg{randomForest}).
#' @param feature_df A data frame containing all required input features, including:
#'   \code{gc_content}, \code{RNA_type}, \code{RNA_region}, \code{exon_length},
#'   \code{distance_to_junction}, \code{evolutionary_conservation}, and \code{DNA_5mer}.
#' @param positive_threshold A numeric cutoff (default = 0.5) to classify predictions as "Positive" or "Negative".
#'
#' @return A data frame that includes all input features plus two additional columns:
#'   \itemize{
#'     \item \code{predicted_m6A_prob}: predicted probability for m6A modification.
#'     \item \code{predicted_m6A_status}: predicted class label ("Positive" or "Negative").
#'   }
#'
#' @examples
#' # prediction_multiple(ml_fit, feature_df)
#'
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  preds <- predict(ml_fit, feature_df, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- preds
  feature_df$predicted_m6A_status <- ifelse(preds > positive_threshold, "Positive", "Negative")
  #Complete this function by writing code at the `___`
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

#' Predict m6A site for a single RNA entry
#'
#' This function predicts m6A modification probability and classification
#' for a single RNA sequence input based on its features.
#'
#' @param ml_fit A trained machine learning model (e.g., from \pkg{randomForest}).
#' @param gc_content Numeric value of GC content for the RNA sequence.
#' @param RNA_type Factor indicating RNA type, e.g., "mRNA", "lincRNA", "lncRNA", or "pseudogene".
#' @param RNA_region Factor indicating RNA region, e.g., "CDS", "intron", "3'UTR", or "5'UTR".
#' @param exon_length Numeric value representing the exon length.
#' @param distance_to_junction Numeric value for the distance to the nearest exon junction.
#' @param evolutionary_conservation Numeric score for evolutionary conservation.
#' @param DNA_5mer A character string representing the 5-mer DNA sequence.
#' @param positive_threshold Numeric value (default = 0.5) used as a cutoff for classification.
#'
#' @return A named vector containing:
#'   \itemize{
#'     \item \code{predicted_m6A_prob}: predicted probability for m6A modification.
#'     \item \code{predicted_m6A_status}: predicted class ("Positive" or "Negative").
#'   }
#'
#' @examples
#' # prediction_single(ml_fit, 0.45, "mRNA", "CDS", 1500, 120, 0.8, "ATCGA")
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  RNA_type <- factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  RNA_region <- factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  feature_df <- data.frame(gc_content = gc_content, RNA_type = RNA_type, RNA_region = RNA_region, exon_length = exon_length, distance_to_junction = distance_to_junction, evolutionary_conservation = evolutionary_conservation, DNA_5mer = DNA_5mer)
  res <- prediction_multiple(ml_fit, feature_df, positive_threshold)
  returned_vector <- c(predicted_m6A_prob = res$predicted_m6A_prob[1], predicted_m6A_status = res$predicted_m6A_status[1])
  #Complete this function by writing code at the `___`
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
