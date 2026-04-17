# Function to calculate top 1% IBS threshold and plot distributions

plot_ibs_distribution <- function(background_data, matched_pairs_data) {
    # Calculate the IBS threshold for top 1%
    threshold_value <- quantile(background_data, 0.99)
    
    # Create a histogram for background DNA-RNA IBS distribution
    hist(background_data, breaks=30, main='Background DNA-RNA IBS Distribution', xlab='IBS Values', col='lightblue', xlim=c(0, max(c(background_data, matched_pairs_data))))
    abline(v=threshold_value, col='red', lwd=2)
    text(x=threshold_value, y=max(hist(background_data, breaks=30, plot=FALSE)$density), label=paste('Top 1% Threshold:', round(threshold_value, 2)), pos=4, col='red')
    
    # Overlay the matched pairs distribution
    hist(matched_pairs_data, breaks=30, add=TRUE, col=rgb(0, 1, 0, 0.5), main='Matched Pairs Distribution')
    
    legend('topright', legend=c('Background Distribution', 'Matched Pairs Distribution', 'Top 1% Threshold'), fill=c('lightblue', rgb(0, 1, 0, 0.5), 'red'))
}

# Example usage:
# plot_ibs_distribution(background_data_vector, matched_pairs_vector)
