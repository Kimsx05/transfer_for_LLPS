function (object, min.cells = 10, min.samples = NULL, rare.keep = FALSE, 
    nonFilter.keep = FALSE) 
{
    net <- object@net
    if (nonFilter.keep == TRUE) {
        cat("The non-filtered cell-cell communication is stored in `object@net$prob.nonFilter` and `object@net$pval.nonFilter`. \n")
        object@net$prob.nonFilter <- net$prob
        object@net$pval.nonFilter <- net$pval
    }
    num.interaction0 <- sum(net$prob > 0)
    cell.excludes <- which(as.numeric(table(object@idents)) <= 
        min.cells)
    if (length(cell.excludes) > 0) {
        cat("The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ", 
            toString(levels(object@idents)[cell.excludes]), "!", 
            "\t")
        net$prob[cell.excludes, , ] <- 0
        net$prob[, cell.excludes, ] <- 0
        num.interaction1 <- sum(net$prob > 0)
        pct.dicrease <- scales::percent((num.interaction0 - num.interaction1)/num.interaction0, 
            accuracy = 0.1)
        cat(paste0(pct.dicrease, " interactions are removed!", 
            "\n"))
    }
    else {
        num.interaction1 <- num.interaction0
    }
    sample.info <- object@meta$samples
    sample.id <- levels(sample.info)
    if (is.null(min.samples)) {
        min.samples <- 1
    }
    else if (min.samples > length(sample.id)) {
        stop(paste0("There are only ", length(sample.id), " samples in the data. Please change the value of `min.samples`! "))
    }
    if (length(sample.id) >= 2 & min.samples >= 2) {
        if (object@options$parameter$raw.use == TRUE) {
            data <- as.matrix(object@data.signaling)
        }
        else {
            if ("data.smooth" %in% methods::slotNames(object) == 
                FALSE) {
                stop("`object@data.smooth` is missing. Please update the CellChat object via `updateCellChat`! \n")
            }
            data <- as.matrix(object@data.smooth)
        }
        data.use <- data/max(data)
        group <- object@idents
        type <- object@options$parameter$type.mean
        trim <- object@options$parameter$trim
        FunMean <- switch(type, triMean = triMean, truncatedMean = function(x) mean(x, 
            trim = trim, na.rm = TRUE), thresholdedMean = function(x) thresholdedMean(x, 
            trim = trim, na.rm = TRUE), median = function(x) median(x, 
            na.rm = TRUE))
        LR <- dimnames(net$prob)[[3]]
        idx.nonzero <- which(apply(net$prob, 3, sum) != 0)
        LR.nonzero <- LR[idx.nonzero]
        interaction_input <- object@DB$interaction
        complex_input <- object@DB$complex
        geneIfo <- object@DB$geneInfo
        idx <- match(LR.nonzero, interaction_input$interaction_name)
        geneL <- as.character(interaction_input$ligand[idx])
        geneR <- as.character(interaction_input$receptor[idx])
        geneLR <- c(unique(geneL), unique(geneR))
        geneLR <- extractGeneSubset(geneLR, complex_input, geneIfo)
        data.use <- data.use[rownames(data.use) %in% geneLR, 
            ]
        score.LR <- array(0, dim = c(nlevels(group), nlevels(group), 
            length(LR.nonzero), length(sample.id)))
        LR.nonzero.all <- c()
        cell.excludes.sample <- c()
        for (i in 1:length(sample.id)) {
            cell.use <- which(sample.info == sample.id[i])
            group.use <- group[cell.use]
            group.use <- droplevels(group.use)
            cell.excludes.sample.i <- which(as.numeric(table(object@idents[cell.use])) <= 
                min.cells)
            cell.excludes.sample <- c(cell.excludes.sample, cell.excludes.sample.i)
            data.use.i <- data.use[, cell.use]
            data.use.avg <- aggregate(t(data.use.i), list(group.use), 
                FUN = FunMean)
            data.use.avg <- t(data.use.avg[, -1])
            group.exist <- which(levels(group) %in% unique(group.use))
            if (length(group.exist) < nlevels(group)) {
                data.use.avg.temp <- matrix(0, nrow = nrow(data.use), 
                  ncol = nlevels(group))
                data.use.avg.temp[, group.exist] <- data.use.avg
                rownames(data.use.avg.temp) <- rownames(data.use.avg)
                data.use.avg <- data.use.avg.temp
            }
            colnames(data.use.avg) <- levels(group)
            dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
            dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
            for (jj in 1:length(LR.nonzero)) {
                score.LR[, , jj, i] <- Matrix::crossprod(matrix(dataLavg[jj, 
                  ], nrow = 1), matrix(dataRavg[jj, ], nrow = 1))
            }
            if (length(cell.excludes.sample.i) > 0) {
                cat(paste0("The number of cells of the following cell groups in ", 
                  sample.id[i], " sample are less than ", min.cells, 
                  " cells: ", toString(levels(object@idents)[cell.excludes.sample.i]), 
                  "!", "\n"))
                score.LR[cell.excludes.sample.i, , , i] <- 0
                score.LR[, cell.excludes.sample.i, , i] <- 0
            }
        }
        cell.excludes.sample <- unique(cell.excludes.sample)
        if (length(cell.excludes.sample) > 0) {
            cell.excludes.sample <- setdiff(cell.excludes.sample, 
                cell.excludes)
        }
        score.LR[score.LR > 0] <- 1
        score.LR.consitent <- array(0, dim = c(nlevels(group), 
            nlevels(group), length(LR.nonzero)))
        LR.inconsitent <- c()
        for (jj in 1:length(LR.nonzero)) {
            score.LR.sum <- apply(score.LR[, , jj, ], c(1, 2), 
                sum)
            if (sum((score.LR.sum > 0) * (score.LR.sum < min.samples)) > 
                0) {
                score.LR.consitent <- (score.LR.sum >= min.samples) * 
                  1
                if (rare.keep == TRUE & length(cell.excludes.sample) > 
                  0) {
                  score.LR.consitent[cell.excludes.sample, ] <- 1
                  score.LR.consitent[, cell.excludes.sample] <- 1
                }
                net$prob[, , LR.nonzero[jj]] <- net$prob[, , 
                  LR.nonzero[jj]] * score.LR.consitent
            }
        }
        num.interaction2 <- sum(net$prob > 0)
        pct.dicrease <- scales::percent((num.interaction1 - num.interaction2)/num.interaction1, 
            accuracy = 0.1)
        cat(paste0(pct.dicrease, " interactions are removed due to their inconsistence across ", 
            min.samples, " samples!", "\n"))
    }
    object@net <- net
    return(object)
}
