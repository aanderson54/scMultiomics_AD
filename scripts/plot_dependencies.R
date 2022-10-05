
Bed_PeakPlot<-function (object, region, group.by = NULL, 
                       color = "dimgrey") 
{
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  if(class(object) =="data.frame"){
  colnames(object)[c(1,2,3)]<-c("chr","start","end")
  peaks<-makeGRangesFromDataFrame(object,
                                keep.extra.columns=T,
                                seqnames.field="chr",
                                start.field="start",
                                end.field="end",
                                strand.field="strand")
  }else{
    peaks<-object
  }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    if (!is.null(x = group.by)) {
      if (!(group.by %in% colnames(x = peak.df))) {
        warning("Requested grouping variable not found")
        group.by <- NULL
      }
      peak.df$start[peak.df$start < start.pos] <- start.pos
      peak.df$end[peak.df$end > end.pos] <- end.pos
      peak.plot <- ggplot(data = peak.df, aes(color = encodeLabel)) + geom_segment(aes(x = start, y = 0,  xend = end, yend = 0), size = 2, data = peak.df)
    }else{
      peak.df$start[peak.df$start < start.pos] <- start.pos
      peak.df$end[peak.df$end > end.pos] <- end.pos
      peak.plot <- ggplot(data = peak.df) + geom_segment(aes(x = start, y = 0, xend = end, yend = 0), size = 2, data = peak.df)
    } 
    
  }
  else {
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() + ylab(label = "Peaks") + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  if (is.null(x = group.by)) {
    peak.plot <- peak.plot + scale_color_manual(values = color) + 
      theme(legend.position = "none")
  }
  return(peak.plot)
}



Bed_GWASPlot<-function (object, region, group.by = NULL, 
                        color = "dimgrey") 
{
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  if(class(object) =="data.frame"){
    colnames(object)[c(1,2,3)]<-c("chr","start","end")
    peaks<-makeGRangesFromDataFrame(object,
                                    keep.extra.columns=T,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end",
                                    strand.field="strand")
  }else{
    peaks<-object
  }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    if (!is.null(x = group.by)) {
      if (!(group.by %in% colnames(x = peak.df))) {
        warning("Requested grouping variable not found")
        group.by <- NULL
      }
      peak.df$start[peak.df$start < start.pos] <- start.pos
      peak.df$end[peak.df$end > end.pos] <- end.pos
      peak.plot <- ggplot(data = peak.df, aes(color = group.by)) + geom_point(aes(x = start, y = 0), size = 3, shape=18,data = peak.df)
    }else{
      peak.df$start[peak.df$start < start.pos] <- start.pos
      peak.df$end[peak.df$end > end.pos] <- end.pos
      peak.plot <- ggplot(data = peak.df) + geom_point(aes(x = start, y = 0), size = 3, shape=18,data = peak.df)
      
    } 
    
  }
  else {
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() + ylab(label = "Peaks") + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  if (is.null(x = group.by)) {
    peak.plot <- peak.plot + scale_color_manual(values = color) + 
      theme(legend.position = "none")
  }
  return(peak.plot)
}





FindRegion <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # if separators are present in the string and we can convert the
    # start to a number, assume we're using genomic coordinates
    if (all(sapply(X = sep, FUN = grepl, x = region))) {
      region <- StringToGRanges(regions = region, sep = sep)
    } else {
      region <- LookupGeneCoords(object = object, assay = assay, gene = region)
      if (is.null(x = region)) {
        stop("Gene not found")
      }
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}





















SingleCoveragePlot <- function(
  object,
  region,
  features = NULL,
  assay = NULL,
  show.bulk = FALSE,
  expression.assay = NULL,
  expression.slot = "data",
  annotation = TRUE,
  peaks = TRUE,
  peaks.group.by = NULL,
  ranges = NULL,
  ranges.group.by = NULL,
  ranges.title = "Ranges",
  region.highlight = NULL,
  links = TRUE,
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  group.by = NULL,
  window = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  heights = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("Requested assay is not a ChromatinAssay.")
  }
  if (!is.null(x = group.by)) {
    Idents(object = object) <- group.by
  }
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  region <- FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  reads.per.group <- AverageCounts(
    object = object,
    group.by = group.by,
    verbose = FALSE
  )
  cells.per.group <- CellsPerGroup(
    object = object,
    group.by = group.by
  )
  cutmat <- CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- start(x = region):end(x = region)
  group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
  scale.factor <- SetIfNull(
    x = scale.factor, y = median(x = group.scale.factors)
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  p <- CoverageTrack(
    cutmat = cutmat,
    region = region,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    window = window,
    ymax = ymax,
    obj.groups = obj.groups,
    region.highlight = region.highlight,
    downsample.rate = downsample.rate,
    max.downsample = max.downsample
  )
  # create bigwig tracks
  if (!is.null(x = bigwig)) {
    if (!inherits(x = bigwig, what = "list")) {
      warning("BigWig should be a list of file paths")
      bigwig <- list("bigWig" = bigwig)
    }
    if (length(x = bigwig.type == 1)) {
      bigwig.type <- rep(x = bigwig.type, length(x = bigwig))
    } else if (length(x = bigwig.type) != length(x = bigwig)) {
      stop("Must supply a bigWig track type for each bigWig file")
    }
    # iterate over list of files
    bigwig.tracks <- list()
    for (i in seq_along(along.with = bigwig)) {
      bigwig.tracks[[i]] <- BigwigTrack(
        region = region,
        bigwig = bigwig[[i]],
        y_label = names(x = bigwig)[[i]],
        type = bigwig.type[[i]]
      )
    }
  } else {
    bigwig.tracks <- NULL
  }
  if (!is.null(x = features)) {
    ex.plot <- ExpressionPlot(
      object = object,
      features = features,
      assay = expression.assay,
      idents = idents,
      group.by = group.by,
      slot = expression.slot
    )
    widths <- c(10, length(x = features))
  } else {
    ex.plot <- NULL
    widths <- NULL
  }
  if (annotation) {
    gene.plot <- AnnotationPlot(object = object[[assay]], region = region)
  } else {
    gene.plot <- NULL
  }
  if (links) {
    link.plot <- LinkPlot(object = object[[assay]], region = region)
  } else {
    link.plot <- NULL
  }
  if (peaks) {
    peak.plot <- PeakPlot(
      object = object,
      assay = assay,
      region = region,
      group.by = peaks.group.by
    )
  } else {
    peak.plot <- NULL
  }
  if (!is.null(x = ranges)) {
    range.plot <- PeakPlot(
      object = object,
      assay = assay,
      region = region,
      peaks = ranges,
      group.by = ranges.group.by,
      color = "brown3") +
      ylab(ranges.title)
  } else {
    range.plot <- NULL
  }
  if (tile) {
    # reuse cut matrix
    tile.df <- ComputeTile(
      cutmatrix = cutmat,
      groups = obj.groups,
      window = tile.size,
      n = tile.cells,
      order = "total"
    )
    tile.plot <- CreateTilePlot(
      df = tile.df,
      n = tile.cells
    )
  } else {
    tile.plot <- NULL
  }
  if (show.bulk) {
    object$bulk <- "All cells"
    reads.per.group <- AverageCounts(
      object = object,
      group.by = "bulk",
      verbose = FALSE
    )
    cells.per.group <- CellsPerGroup(
      object = object,
      group.by = "bulk"
    )
    bulk.scale.factor <- suppressWarnings(reads.per.group * cells.per.group)
    bulk.groups <- rep(x = "All cells", length(x = obj.groups))
    names(x = bulk.groups) <- names(x = obj.groups)
    
    bulk.plot <- CoverageTrack(
      cutmat = cutmat,
      region = region,
      group.scale.factors = bulk.scale.factor,
      scale.factor = scale.factor,
      window = window,
      ymax = ymax,
      obj.groups = bulk.groups,
      downsample.rate = downsample.rate,
      max.downsample = max.downsample
    ) +
      scale_fill_manual(values = "grey") +
      ylab("")
  } else {
    bulk.plot <- NULL
  }
  nident <- length(x = unique(x = obj.groups))
  bulk.height <- (1 / nident) * 10
  bw.height <- bulk.height * length(x = bigwig.tracks)
  heights <- SetIfNull(
    x = heights, y = c(10, bulk.height, bw.height, 10, 2, 1, 1, 3)
  )
  # pre-combine bigwig tracks
  if (!is.null(x = bigwig.tracks)) {
    bigwig.tracks <- CombineTracks(
      plotlist = bigwig.tracks,
      heights = rep(x = 1, length(x = bigwig.tracks))
    )
  }
  p <- CombineTracks(
    plotlist = list(p, bulk.plot, bigwig.tracks, tile.plot, gene.plot,
                    peak.plot, range.plot, link.plot),
    expression.plot = ex.plot,
    heights = heights,
    widths = widths
  ) & theme(
    legend.key.size = unit(x = 1/2, units = "lines"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
  return(p)
}








CoverageTrack <- function(
  cutmat,
  region,
  group.scale.factors,
  scale.factor,
  obj.groups,
  ymax,
  downsample.rate,
  region.highlight = NULL,
  window = 100,
  max.downsample = 3000
) {
  window.size <- width(x = region)
  levels.use <- levels(x = obj.groups)
  coverages <- ApplyMatrixByGroup(
    mat = cutmat,
    fun = colSums,
    groups = obj.groups,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    normalize = TRUE
  )
  if (!is.na(x = window)) {
    coverages <- group_by(.data = coverages, group)
    coverages <- mutate(.data = coverages, coverage = roll_sum(
      x = norm.value, n = window, fill = NA, align = "center"
    ))
    coverages <- ungroup(x = coverages)
  } else {
    coverages$coverage <- coverages$norm.value
  }
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  coverages <- coverages[!is.na(x = coverages$coverage), ]
  coverages <- group_by(.data = coverages, group)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = coverages, n = sampling)
  
  # restore factor levels
  if (!is.null(x = levels.use)) {
    colors_all <- hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
    coverages$group <- factor(x = coverages$group, levels = levels.use)
  }
  ymax <- SetIfNull(x = ymax, y = signif(
    x = max(coverages$coverage, na.rm = TRUE), digits = 2)
  )
  ymin <- 0
  
  gr <- GRanges(
    seqnames = chromosome,
    IRanges(start = start.pos, end = end.pos)
  )
  p <- ggplot(
    data = coverages,
    mapping = aes(x = position, y = coverage, fill = group)
  ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("Normalized accessibility \n(range ",
                        as.character(x = ymin), " - ",
                        as.character(x = ymax), ")")) +
    ylim(c(ymin, ymax)) +
    theme_browser(legend = FALSE) +
    theme(panel.spacing.y = unit(x = 0, units = "line"))
  if (!is.null(x = levels.use)) {
    p <- p + scale_fill_manual(values = colors_all)
  }
  if (!is.null(x = region.highlight)) {
    if (!inherits(x = region.highlight, what = "GRanges")) {
      warning("region.highlight must be a GRanges object")
    } else {
      md <- mcols(x = region.highlight)
      if ("color" %in% colnames(x = md)) {
        color.use <- md$color
      } else {
        color.use <- rep(x = "grey", length(x = region.highlight))
      }
      df <- data.frame(
        "start" = start(x = region.highlight),
        "end" = end(x = region.highlight),
        "color" = color.use
      )
      df$start <- ifelse(
        test = df$start < start.pos,
        yes = start.pos,
        no = df$start
      )
      df$end <- ifelse(
        test = df$end > end.pos,
        yes = end.pos,
        no = df$end
      )
      p <- p +
        geom_rect(
          data = df,
          inherit.aes = FALSE,
          aes_string(
            xmin = "start",
            xmax = "end",
            ymin = 0,
            ymax = ymax),
          fill = rep(x = df$color, length(x = unique(x = coverages$group))),
          color = "transparent",
          alpha = 0.2
        )
    }
  }
  return(p)
}





BigwigTrack <- function(
  region,
  bigwig,
  smooth = 200,
  type = "coverage",
  y_label = "Score",
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  possible_types <- c("line", "heatmap", "coverage")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  if (.Platform$OS.type == "windows") {
    message("BigwigTrack not supported on Windows")
    return(NULL)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  if (!requireNamespace("RcppRoll", quietly = TRUE)) {
    message("Please install RcppRoll. install.packages('RcppRoll')")
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  region_data <- rtracklayer::import(
    con = bigwig,
    which = region,
    as = "NumericList"
  )[[1]]
  if (!is.null(x = smooth)) {
    region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
  }
  region_data <- data.frame(
    position = start(x = region):end(x = region),
    score = region_data,
    stringsAsFactors = FALSE
  )
  window.size = width(x = region)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = region_data, n = sampling)
  p <- ggplot(
    data = coverages,
    mapping = aes_string(x = "position", y = "score")
  )
  if (type == "line") {
    p <- p + geom_line()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    region_data$bin <- floor(x = region_data$position / smooth)
    region_data <- group_by(region_data, bin)
    region_data <- mutate(region_data, score = mean(x = score))
    region_data <- ungroup(region_data)
    region_data <- unique(x = region_data[, c("bin", "score")])
    p <- ggplot(
      data = region_data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c()
  } else if (type == "coverage") {
    p <- p + geom_area(aes(fill="A"))+theme(legend.position="none")
  }
  chromosome <- as.character(x = seqnames(x = region))
  p <- p + theme_classic() +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label)
  return(p)
}



LinkPlot.height <- function(object, region, min.cutoff = -1.5, max.cutoff=1.7) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  chromosome <- seqnames(x = region)
  
  # extract link information
  links <- Links(object = object)
  
  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }
  
  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)
  
  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  
  # remove links outside region
  link.df <- link.df[link.df$start >= start(x = region) & link.df$end <= end(x = region), ]
  
  # plot
  if (nrow(x = link.df) > 0) {
    # convert to format for geom_bezier
    link.df$group <- seq_len(length.out = nrow(x = link.df))
    df <- data.frame(
      x = c(link.df$start,
            (link.df$start + link.df$end) / 2,
            link.df$end),
      y = c(rep(x = 0, nrow(x = link.df)),
            rep(x = -1, nrow(x = link.df)),
            rep(x = 0, nrow(x = link.df))),
      group = rep(x = link.df$group, 3),
      score = rep(link.df$group2, 3),
      height=rep(link.df$score*2,3)
    )
    df$y<-ifelse(df$y==-1, df$height, df$y)
    tmp<-data.frame(x = c(link.df$start[1],
                          (link.df$start[1] + link.df$end[1]) / 2,
                          link.df$end[1]),
                    y = c(0,max.cutoff,0),
                    group = rep(max(link.df$group)+1, 3),
                    score = c(10,10,10),
                    height=c(0,0,0))
    if (nrow(df[which(df$y<0),])>0){
      tmp2<-data.frame(x = c(link.df$start[1],
                             (link.df$start[1] + link.df$end[1]) / 2,
                             link.df$end[1]),
                       y = c(0,min.cutoff,0),
                       group = rep(max(link.df$group)+2, 3),
                       score = c(10,10,10),
                       height=c(0,0,0))
      tmp<-rbind(tmp,tmp2)
    }
    df<-rbind(df, tmp)
    p <- ggplot(data = df) +
      geom_bezier(
        mapping = aes_string(x = "x", y = "y", group = "group", color = "score")
      ) +scale_color_gradient2(low="green",mid="orange",high="grey40", midpoint=1, na.value=alpha("white", 0),limits=c(0, 2))+
      geom_hline(yintercept = 0, color = 'grey') 
  } else {
    p <- ggplot(data = link.df)
  }
  p <- p +
    theme_classic() +
    theme(legend.position="none") +
    ylab("Links") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
  return(p)
}


LinkPlot.height.Perm <- function(object, region, min.cutoff = -1.5, max.cutoff=1.7) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  chromosome <- seqnames(x = region)
  
  # extract link information
  links <- Links(object = object)
  
  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }
  
  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)
  
  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  
  # remove links outside region
  link.df <- link.df[link.df$start >= start(x = region) & link.df$end <= end(x = region), ]
  
  # plot
  if (nrow(x = link.df) > 0) {
    # convert to format for geom_bezier
    link.df$group <- seq_len(length.out = nrow(x = link.df))
    df <- data.frame(
      x = c(link.df$start,
            (link.df$start + link.df$end) / 2,
            link.df$end),
      y = c(rep(x = 0, nrow(x = link.df)),
            rep(x = -1, nrow(x = link.df)),
            rep(x = 0, nrow(x = link.df))),
      group = rep(x = link.df$group, 3),
      score = rep(link.df$group2, 3),
      height=rep(link.df$score*2,3)
    )
    df$y<-ifelse(df$y==-1, df$height, df$y)
    df$alpha=1
    tmp<-data.frame(x = c(link.df$start[1],
                          (link.df$start[1] + link.df$end[1]) / 2,
                          link.df$end[1]),
                    y = c(0,max.cutoff,0),
                    group = rep(max(link.df$group)+1, 3),
                    score = c(1000,1000,1000),
                    height=c(0,0,0))
      tmp2<-data.frame(x = c(link.df$start[1],
                             (link.df$start[1] + link.df$end[1]) / 2,
                             link.df$end[1]),
                       y = c(0,min.cutoff,0),
                       group = rep(max(link.df$group)+2, 3),
                       score = c(1000,1000,1000),
                       height=c(0,0,0))
      tmp<-rbind(tmp,tmp2)
      tmp$alpha=0
    
    df<-rbind(df, tmp)
    df$alpha<-as.factor(df$alpha)
    p <- ggplot(data = df) +
      geom_bezier(
        mapping = aes_string(x = "x", y = "y", group = "group", color = "score", alpha="alpha")
      ) +scale_color_viridis(limits=c(0,100), direction=-1, na.value="white")+
      geom_hline(yintercept = 0, color = 'grey') +scale_alpha_discrete(range=c(0,1))
  } else {
    p <- ggplot(data = link.df)
  }
  p <- p +
    theme_classic() + 
    ylab("Links") +
    xlab(label = paste0(chromosome, "position (bp)")) +
    xlim(c(start(x = region), end(x = region)))+guides(alpha = "none")
  return(p)
}



LinkPlot.height.alpha <- function(object, region, min.cutoff = -1.5, max.cutoff=1.7) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  chromosome <- seqnames(x = region)
  
  # extract link information
  links <- Links(object = object)
  
  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }
  
  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)
  
  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  
  # remove links outside region
  link.df <- link.df[link.df$start >= start(x = region) & link.df$end <= end(x = region), ]
  
  # plot
  if (nrow(x = link.df) > 0) {
    # convert to format for geom_bezier
    link.df$group <- seq_len(length.out = nrow(x = link.df))
    df <- data.frame(
      x = c(link.df$start,
            (link.df$start + link.df$end) / 2,
            link.df$end),
      y = c(rep(x = 0, nrow(x = link.df)),
            rep(x = -1, nrow(x = link.df)),
            rep(x = 0, nrow(x = link.df))),
      group = rep(x = link.df$group, 3),
      score = rep(link.df$group2, 3),
      height=rep(link.df$score*2,3),
      alpha=rep(link.df$alpha,3)
    )
    df$y<-ifelse(df$y==-1, df$height, df$y)
    tmp<-data.frame(x = c(link.df$start[1],
                          (link.df$start[1] + link.df$end[1]) / 2,
                          link.df$end[1]),
                    y = c(0,max.cutoff,0),
                    group = rep(max(link.df$group)+1, 3),
                    score = c(10,10,10),
                    height=c(0,0,0),
                    alpha=c(0,0,0))
    if (nrow(df[which(df$y<0),])>0){
      tmp2<-data.frame(x = c(link.df$start[1],
                             (link.df$start[1] + link.df$end[1]) / 2,
                             link.df$end[1]),
                       y = c(0,min.cutoff,0),
                       group = rep(max(link.df$group)+2, 3),
                       score = c(10,10,10),
                       height=c(0,0,0),
                       alpha=c(0,0,0))
      tmp<-rbind(tmp,tmp2)
    }
    df<-rbind(df, tmp)
    p <- ggplot(data = df) +
      geom_bezier(
        mapping = aes_string(x = "x", y = "y", group = "group", color = "score", alpha="alpha")
      ) +scale_color_gradient2(low="red",mid="grey54",high="blue", midpoint=1, na.value=alpha("white", 0),limits=c(0, 2))+
      geom_hline(yintercept = 0, color = 'grey') 
  } else {
    p <- ggplot(data = link.df)
  }
  p <- p +
    theme_classic() +
    theme(
      legend.position="none") +
    ylab("Links") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
  return(p)
}









AnnotationPlot2 <- function(object, region, mode = "gene") {
  if(mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 3
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.2
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = 1.5, angle=45
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}








