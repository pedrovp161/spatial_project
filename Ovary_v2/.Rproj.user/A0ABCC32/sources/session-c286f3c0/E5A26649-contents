
library(tidyverse)
library(Seurat)

# functions

get_clu_barcode <- function(sp, sample, clu) {
    
    bc <-  sp@meta.data %>%
        rownames_to_column(var = 'barcode') %>%
        as_tibble() %>%
        filter(Sample == sample) %>%
        filter(seurat_clusters == clu) %>%
        pull(barcode)
    
    return(bc)
}

get_clu_barcode(ov_combined, 'ER_1', 0)
GetTissueCoordinates(ov_combined, "ER_1")



clu_dist <- function(sp, sample, img, clu_1, clu_2) {
    
    coord_1 <- GetTissueCoordinates(sp, img)[get_clu_barcode(sp, sample, clu_1),]
    coord_2 <- GetTissueCoordinates(sp, img)[get_clu_barcode(sp, sample, clu_2),]
    
    d <- c()
    
    for (i in 1:nrow(coord_1)) {
        d[i] <- median(sqrt((coord_1[i,1] - coord_2[,1])^2 + (coord_1[i,2] - coord_2[,2])^2))
    }
    
    return(d)

}

clu_dist(ov_combined, sample = 'ER_1', img = 'ER_1', clu_1 = 0, clu_2 = 1)




dist <- list()
len <- 0

for (s in ov_combined@meta.data$Sample %>% unique()) {
    for (i in 0:8) {
        for (j in 0:8) {
            len <- len + 1
            dist[[len]] <- tibble(sample = s,
                                  cluster_1 = as.character(i),
                                  cluster_2 = as.character(j),
                                  distance = clu_dist(ov_combined, sample = s, img = s, clu_1 = i, clu_2 = j))
        }
    }
}

dist <- do.call("rbind", dist)

dist <- dist %>%
    mutate(group = ifelse(str_detect(sample, 'ER'), 'ER', 'PR')) %>%
    group_by(group, cluster_1, cluster_2) %>%
    summarise(m_distance = median(distance, na.rm = T))


# symmetry


dist <- dist %>%
    mutate(Key = ifelse(cluster_1 < cluster_2, 
                        paste0(cluster_1, '_', cluster_2),
                        paste0(cluster_2, '_', cluster_1)))


pdf('./figures/Spot_Distance_052621.pdf', width = 13, height = 7)
dist %>%
    group_by(group, Key) %>%
    summarise(dist = mean(m_distance)) %>%
    mutate(Cluster_1 = str_extract(Key, "^[0-9]"),
           Cluster_2 = str_extract(Key, "[0-9]$")) %>%
    filter(group == 'ER') %>%
    drop_na() %>%
    ggplot() +
    geom_tile(aes(x = Cluster_1, y = Cluster_2, fill = dist)) +
    coord_equal() +
    theme_bw() +
    ggtitle('ER') +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
dist %>%
    group_by(group, Key) %>%
    summarise(dist = mean(m_distance)) %>%
    mutate(Cluster_1 = str_extract(Key, "^[0-9]"),
           Cluster_2 = str_extract(Key, "[0-9]$")) %>%
    filter(group == 'PR') %>%
    drop_na() %>%
    ggplot() +
    geom_tile(aes(x = Cluster_1, y = Cluster_2, fill = dist)) +
    coord_equal() +
    theme_bw() +
    ggtitle('PR') +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    plot_layout(guides = "collect") &
    scale_fill_distiller(palette = 'RdYlBu', direction = 1, limits = range(c(0,120)))
dev.off()

dist %>%
    group_by(group, Key) %>%
    summarise(dist = mean(m_distance)) %>%
    mutate(Cluster_1 = str_extract(Key, "^[0-9]"),
           Cluster_2 = str_extract(Key, "[0-9]$")) %>%
    filter(group == 'PR') %>%
    drop_na() %>%
    print(n = 34)
    


# 032421

pdf('./figures/Spot_Distance_042121.pdf', width = 13, height = 7)
dist %>%
    filter(group == 'ER') %>%
    drop_na() %>%
    ggplot() +
    geom_tile(aes(x = cluster_1, y = cluster_2, fill = m_distance)) +
    scale_fill_distiller(palette = 'RdYlBu', direction = 1) +
    coord_equal() +
    theme_bw() +
    ggtitle('ER') +
dist %>%
    filter(group == 'PR') %>%
    filter(cluster_1 != '6' & cluster_2 != '6') %>%
    ggplot() +
    geom_tile(aes(x = cluster_1, y = cluster_2, fill = m_distance)) +
    scale_fill_distiller(palette = 'RdYlBu', direction = 1) +
    coord_equal() +
    theme_bw() +
    ggtitle('PR')
dev.off()


for (i in 0:6) {
    
    dist[[i+1]] <- tibble(cluster = as.character(i),
                        distance = clu_dist(ov, sample = 'Poor', img = 'P1', clu_1 = 1, clu_2 = i))
    
}

dist <- do.call("rbind", dist)

dist %>%
    filter(cluster != '1') %>%
    ggplot(aes(x = cluster, y = distance)) +
    geom_boxplot(aes(fill = cluster), alpha = 0.5) +
    geom_jitter(alpha = 0.2, width = 0.15) +
    scale_fill_brewer(palette = "Set2") +
    ggtitle("Median Distance with Cluster 1")



