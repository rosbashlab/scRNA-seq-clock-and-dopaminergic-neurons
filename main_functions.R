
# save figures
save_pdf <- function(x){
  paste0("Figure/",x,"_",Sys.Date(),".pdf")
}


# calculate entropy in each cell

f.entropy <- function( m.in ){
  m.in <- apply( m.in, 2, function( x ){ x / sum(x, na.rm=T ); } );
  v.entropy <- apply( m.in, 2, function(x){ y <- unname(x); y <- y[ y > 0 ]; sum( -y * log(y) ); } );
  return( v.entropy );
}


# integration and clustering

f.cluster.seuratv3.2 <- function( s.in, nDims = 50, k.anchor.use = 5, k.filter.use = 199, k.score.use = 30, k.weight.use = 100, split.by.var = "time" ){
  
  cat( "f.cluster.seurat.v3:\n" );
  cat( "k.anchor.use=", k.anchor.use, "k.filter.use=", k.filter.use, "k.score.use=", k.score.use, "k.weight.use=", k.weight.use, "split.by.var=", split.by.var, "\n" );
  
  
  v.mt.genes <- grep( "^mt:", rownames( s.in@assays$RNA@data ), value=T, ignore.case=T );
  v.ribo.genes <- grep( "^rp[ls]", rownames( s.in@assays$RNA@data ), value=T, ignore.case=T );
  v.rRNA.genes <- grep( "rRNA", rownames( s.in@assays$RNA@data ), value=T, ignore.case=T );
  v.tRNA.genes <- grep( "tRNA", rownames( s.in@assays$RNA@data ), value=T, ignore.case=T );
  v.ercc.genes <- grep( "^ERCC", rownames( s.in@assays$RNA@data ), value=T, ignore.case=F);
  v.var.genes.exclude <- mixedsort( unique( c( v.mt.genes, v.ribo.genes, v.rRNA.genes, v.tRNA.genes, v.ercc.genes, "GFP" ) ) );
  
  remove(v.mt.genes, v.ribo.genes, v.rRNA.genes,v.tRNA.genes,v.ercc.genes)
  
  cat( "there are ", length( v.var.genes.exclude)," genes will be excluded for clustering\n\n")
  
  cat( "splitting data by ", split.by.var, "\n" );
  l.in <- SplitObject( s.in, split.by = split.by.var );
  cat( "split into", length( l.in ), "objects\n" );
  
  cat( "running SCTransform\n" );
  
  for (i in 1:length( l.in ) ) {
    
    l.in[[i]] <- SCTransform( l.in[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"), assay="RNA", verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 3000 );
    
    l.in[[i]]@assays$SCT@var.features <- setdiff( l.in[[i]]@assays$SCT@var.features, v.var.genes.exclude );
    
  }
  
  v.common <- rownames( l.in[[1]]$SCT@data );
  for (  i in 2:length( l.in ) ){
    v.common <- intersect( v.common, rownames( l.in[[i]]$SCT@data ) );
  } 
  length( v.common );
  cat( "length(v.common)=", length(v.common), "\n" );
  
  v.inte.features <- vector();
  if ( 1 ){                                     # find variable features/genes common to all
    v.var.common <- l.in[[1]]@assays$SCT@var.features;
    for (  i in 2:length( l.in ) ){
      v.var.common <- c( v.var.common, l.in[[i]]@assays$SCT@var.features );
    }
    length( unique( v.var.common ) );
    v.var.common <- names( which( table( v.var.common ) == length( l.in ) ) );
    
    
    v.var.common <- setdiff( v.var.common, v.var.genes.exclude );
    length( v.var.common );
    cat( "length(v.var.common)=", length(v.var.common), "\n" );
    v.inte.features <- v.var.common;
  } else {
    v.inte.features <- SelectIntegrationFeatures( object.list = l.in, nfeatures = 3000 );
  }
  
  cat( "length( v.inte.features ) = ", length( v.inte.features ), "\n" );
  
  print( date() );
  cat( "running FindIntegrationAnchors\n" );
  cat( "k.anchor.use=", k.anchor.use, "k.filter.use=", k.filter.use, "k.score.use=", k.score.use, "\n" );
  l.in.anchors <- FindIntegrationAnchors( object.list = l.in, dims = 1:nDims, assay=rep( "SCT", length( l.in )), anchor.features=v.inte.features, k.anchor=k.anchor.use, k.filter=k.filter.use, k.score=k.score.use, verbose=T );
  
  cat( "integrating data:\n\n" );
  cat( "k.weight.use=", k.weight.use, "\n" );
  
  cat( "running IntegrateData\n" );
  s.in.inte <- IntegrateData( anchorset = l.in.anchors, dims = 1:nDims, features.to.integrate=v.common, k.weight = k.weight.use, verbose = T );
  cat( "Data integrated\n" );
  print( date() );
  VariableFeatures( s.in.inte, assay="integrated" ) <- v.inte.features;
  cat( "length(s.in.inte$integrated@var.features)=", length(s.in.inte$integrated@var.features), "\n" );
  
  DefaultAssay( s.in.inte ) <- "integrated";
  s.in.inte <- ScaleData( s.in.inte, verbose=TRUE, assay="integrated", features=rownames( s.in.inte$integrated@data ) );
  cat( "Data scaled\n" );
  
  
  cat( "Running PCA\n" );
  s.in.inte <- RunPCA( s.in.inte, assay="integrated", npcs = nDims, verbose = FALSE, seed.use = 123); cat( "RunPCA done\n" );
  cat( "Running RunTSNE\n" );
  s.in.inte <- RunTSNE( s.in.inte, reduction = "pca", dims = 1:nDims, seed.use = 123 ); cat( "RunTSNE done\n" );
  
  cat( "PCA and TSNE ran\n" );
  
  cat( "Finding Neighbors\n" );
  s.in.inte <- FindNeighbors(s.in.inte, assay="integrated", reduction = "pca", dims = 1:nDims, force.recalc=T );
  cat( "Finding Clusters\n" );
  s.in.inte <- FindClusters(s.in.inte, assay="integrated", resolution = 1.0 )
  
  cat( "table of clusters:", table( s.in.inte@active.ident ), "\n" );
  cat( "table of time: ", table( s.in.inte@meta.data$time ), "\n" );
  cat( "done with f.cluster.seuratv3!\n" );
  return( s.in.inte );
}


f.plot.marker.heatmap <- function(s.in = s.clk.sub,df.marker = df.m.clk,features=v.intersect.tfs){
  
  require(circlize)
  col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "red"))
  
  df.loc <- df.marker[ df.marker$p_val_adj < 0.05 & df.marker$avg_logFC > log(1.25) & df.marker$gene %in% features , ];
  tmp <- sapply( mixedsort(unique(df.marker$cluster)), function(x){ v <- df.loc$gene[ df.loc$cluster == x ]; v <- v[ 1:min(5,length(v) )]; return(v); } )
  str( tmp );
  
  v.features <- unique( unlist( as.vector( tmp ) ) ) %>% as.character();
  
  s.in <- ScaleData(s.in, features = rownames(s.in@assays$RNA@data),assay = "RNA")
  
  v.features <- v.features[v.features %in% rownames(s.clk.sub@assays$RNA@scale.data)]
  
  s.in$seurat_clusters <- factor(s.in$seurat_clusters, levels = mixedsort(unique(as.character(s.in$seurat_clusters))))
  
  df <- sapply(split(s.in$cell.names,s.in$seurat_clusters),function(cells){
    apply(s.in@assays$RNA@scale.data[v.features,cells],1,mean)
  })
  
  df <- df[,mixedsort(colnames(df))];
  
  
  p <- ComplexHeatmap::Heatmap(df,cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               name = "Expression",
                               row_names_gp = gpar(fontsize = 8,fontface = "italic"),
                               border = T,
                               col = col_fun,
                               )
  
  return( p )
  
}


f.feature.tsne <- function(Gene, s.in,reduction = "tsne"){
  
  if( !(Gene %in% rownames(s.in@assays$RNA@data))){
    cat( Gene, "is not identified, double check the gene name!", "\n")
  }
  
  if (reduction == "tsne"){
    df <- as.data.frame( s.in@reductions$tsne@cell.embeddings );
    df$TP10K <- exp( s.in$RNA@data[ Gene, ] ) -1;
    pf <- ggplot( df, aes(x=tSNE_1, y=tSNE_2));
    pf <- pf + geom_point( data=df, aes( x=tSNE_1, y=tSNE_2 ), colour="grey75", size=1.0, alpha=0.5, shape=16 )
    
    pf <- pf + geom_point( data=df[ df$TP10K > 0, ], aes( x=tSNE_1, y=tSNE_2, color=TP10K ), size=1.0, alpha=0.5, shape=16 ) +
      scale_color_gradientn( colors=c("grey75", "red1", "red2", "red3", "red4", "black" ) );
    pf <- pf + labs( title=Gene ) + guides(fill=guide_legend(title="TP10K")) + theme_cowplot();
  } else if (reduction == "umap") {
    df <- as.data.frame( s.in@reductions$umap@cell.embeddings );
    df$TP10K <- exp( s.in$RNA@data[ Gene, ] ) -1;
    pf <- ggplot( df, aes(x=UMAP_1, y=UMAP_2));
    pf <- pf + geom_point( data=df, aes( x=UMAP_1, y=UMAP_2 ), colour="grey75", size=1.0, alpha=0.5, shape=16 )
    
    pf <- pf + geom_point( data=df[ df$TP10K > 0, ], aes( x=UMAP_1, y=UMAP_2, color=TP10K ), size=1.0, alpha=0.5, shape=16 ) +
      scale_color_gradientn( colors=c("grey75", "red1", "red2", "red3", "red4", "black" ) );
    pf <- pf + labs( title=Gene ) + guides(fill=guide_legend(title="TP10K")) + theme_cowplot();
  }
  
  return( pf )
}


f.feature.tsne.v2 <- function(Gene, s.in, threshold = 1) {
  
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(cowplot)
  require(ggrepel)
  require(ggnewscale)
  
  
  if( all(Gene %in% rownames(s.in@assays$RNA@data)) == FALSE ){
    cat( Gene, "is not identified, double check the gene name!", "\n")
  }
  
  
  if( DefaultAssay(s.in) != "RNA" ){    
    DefaultAssay(s.in) <- "RNA"    
  }
  
  l.genes <- Gene;
  
  cat("Genes to plot: ", l.genes, "\n")
  
  df <- as.data.frame( s.in@reductions$tsne@cell.embeddings );
  
  df$ident=(Idents(s.in)[])
  
  df.tp10k <- expm1(FetchData(s.all.subset,vars = l.genes,slot = "data")) %>% 
    set_colnames(c("one","two","three")[1:length(l.genes)])
  
  df <- cbind( df, df.tp10k[rownames(df),])
  
  df_mini <- df[!duplicated(df$ident),]
  
  
  if( length(l.genes) == 1 ){
    
    p <- ggplot(mapping = aes(tSNE_1, tSNE_2)) +
      geom_point( data=df, colour="grey75", size=1.0, alpha=0.5, shape=16 )+
      geom_point( data=df[ df$one > threshold, ], aes(color= one) , size=1.0, alpha=0.5, shape=16 ) + 
      scale_color_gradientn( colors=c("grey75", "red","red1", "red2", "red3", "red4", "black" ), name = l.genes[1])+theme_cowplot()
    
    
  } else if ( length(l.genes) == 2 ) {
    
    p <- ggplot(mapping = aes(tSNE_1, tSNE_2)) +
      geom_point( data=df, colour="grey75", size=1.0, alpha=0.5, shape=16 )+
      geom_point( data=df[ df$one > threshold, ], aes(color= one) , size=1.0, alpha=0.5, shape=16 ) + 
      scale_color_gradientn( colors=c("grey75", "red","red1", "red2", "red3", "red4", "black" ), name = l.genes[1]) +
      new_scale_color()+
      geom_point( data=df[ df$two > threshold, ], aes(color= two), size=1.0, alpha=0.5, shape=16 ) + 
      scale_color_gradientn( colors=c("grey75", "purple", "purple1", "purple2", "purple3", "purple4","Black"), name = l.genes[2])+ theme_cowplot()
    
    
  } else if (length(l.genes) == 3) {
    
    p <-  ggplot(mapping = aes(tSNE_1, tSNE_2)) +
      geom_point( data=df, colour="grey75", size=1.0, alpha=0.5, shape=16 )+
      geom_point( data=df[ df$one > threshold, ], aes(color= one) , size=1.0, alpha=0.5, shape=16 ) + 
      scale_color_gradientn( colors=c("grey75", "red","red1", "red2", "red3", "red4", "black" ), name = l.genes[1]) +
      new_scale_color()+
      geom_point( data=df[ df$two > threshold, ], aes(color= two), size=1.0, alpha=0.5, shape=16 ) + 
      scale_color_gradientn( colors=c("grey75", "purple", "purple1", "purple2", "purple3", "purple4","Black"), name = l.genes[2])+
      new_scale_color()+
      geom_point( data=df[ df$three > threshold, ], aes(color= three), size=1.0, alpha=0.5, shape=16 ) + 
      scale_color_gradientn( colors=c("grey75", "orange", "orange1", "orange2", "orange3", "orange4","Black"), name = l.genes[3])+ theme_cowplot()   
  }
  
  cat("Mission complete \n")
  
  return( p )
  
  
}


# custom heatmap

f.dotplot.custom <- function( p ) {
  
  df<- p$data
  
  exp_mat<-df %>% 
    select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame() 
  
  row.names(exp_mat) <- exp_mat$features.plot  
  exp_mat <- exp_mat[,-1] %>% as.matrix()
  
  percent_mat<-df %>% 
    select(-avg.exp, -avg.exp.scaled) %>%  
    pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame() 
  
  row.names(percent_mat) <- percent_mat$features.plot  
  percent_mat <- percent_mat[,-1] %>% as.matrix()
  
  require(viridis)
  require(Polychrome)

  col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])
  
  cell_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
                gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
  
  
  p1 <- Heatmap(exp_mat,
                heatmap_legend_param=list(title="expression"),
                column_title = "clustered dotplot", 
                col=col_fun,
                rect_gp = gpar(type = "none"),
                cell_fun = cell_fun,
                cluster_columns = F,
                cluster_rows = F,
                row_names_gp = gpar(fontsize = 8),
                border = "black")
  
  return( p1 )
  
}
