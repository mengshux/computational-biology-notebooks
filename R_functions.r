options("stringsAsFactors"=FALSE)
options("repos" = c(CRAN = "http://cran.stat.ucla.edu/"))
options('menu.graphics'=FALSE) # do not use tcl interface for library chooser

if (exists(".env")) rm(.env)
.env <- new.env()
.env$"%nin%" <- function(x,y) !(x %in% y)
.env$"%nlike%" <- function(x,y) !(x %like% y)
.env$pst <- paste0
.env$bioconductor <- function(lib) { source("http://www.bioconductor.org/biocLite.R"); biocLite(lib, lib=.libPaths()[2]) }
.env$loadLibs <- function() { for (theLib in c('data.table','tidyverse','broom','magrittr')) library(theLib, character.only=TRUE) }
.env$install.pack <- function(x) install.packages(x, lib=.libPaths()[2], repos="http://cran.cnr.berkeley.edu/", dependencies=TRUE)
.env$RM <- function() {rm(list=ls())}
.env$rm.mine <- function() {detach(.env); rm(.env)}
.env$ls.mine <- function() ls(envir=.env)
.env$LS <- function() list(GLOBAL=ls(envir=.GlobalEnv), PRIVATE=.env$ls.mine())
.env$Q <- function() q("no")
.env$hist20 <- function() history(20)
.env$h <- utils::head
.env$hh <- function(d,w=5,l=5) {
	if (any(class(d) %in% c("matrix","data.frame"))) {
		a <- min(l, nrow(d))
		b <- min(w, ncol(d))
		if (any(class(d) %in% "data.table")) {
			d[1:a,1:b,with=FALSE]
		} else {
			d[1:a,1:b]
		}
	} else {
		n <- min(l, length(d))
		d[1:n]
	}
}
.env$pn <- function(x, n=50) print(x, n=n)
.env$p100 <- function(x) print(x, n=100)
.env$px <- function(x) print(x,n=nrow(x))

.env$epslog <- function(x, base=10, epsilon=0.1) log(x+epsilon, base=base)
.env$not.na <- function(x,y=NA) {
	if (is.na(y)) {
		x <- x[!is.na(x)]
	} else {
		x[is.na(x)] <- y
	}
	return(x)
}
.env$fix.na <- function(x, fix=0.0) { x[is.na(x)] <- fix; x}
.env$tt <- function(x) table(table(x))
.env$which.na <- function(x) which(is.na(x))
.env$which.nna <- function(x) which(!is.na(x))
.env$closest <- function(x,y) {
	distances <- abs(x - y)
	return(which.min(distances))
}
.env$read.cb <- function(sep="\t", header=TRUE, ...) read.table(pipe("pbpaste"),sep=sep,header=header,...)
.env$tsv.cb <- function(...) read_tsv(pipe("pbpaste"),...)

# substitute y with _v_alue where x is in _k_ey
.env$table_map <- function(x, y, k, v) {
	replace <- v[ match(x,k) ]
	replace[ is.na(replace) ] <- y[ is.na(replace) ]
	return(replace)
}
# replace x with mapped value in lookup.table if it exists
.env$lookup <- function(x, lookup.table) ifelse(x %in% names(lookup.table), lookup.table[x], x)

#.env$o <- function(...) system("open .")
.env$set.width <- function(x=200) { options(width=x); cat("window width =", unlist(options('width')), "\n") }
.env$auto.width <- function() { w <- as.integer(system("/usr/bin/tput cols", intern=TRUE, ignore.stderr=TRUE)); if (length(w) > 0 ) .env$set.width(w) }

.env$use.match <- function(a,b) {
	m <- match(a,b)
	out <-	list(
				match = m,
				a.order = which(!is.na(m)),
				b.order = not.na(m)
			)
	if (identical(a,b)) {
		cat("'a' and 'b' are identical")
	} else if (all(diff(out$a.order) == 1)) {
		cat("All A matched in order\n")
	}
	return(out)
}

.env$gt <- function(threshold, values=values) sum(values > threshold) / length(values)
.env$lt <- function(threshold, values=values) sum(values < threshold) / length(values)

.env$geo.mean <- function(x) exp(mean(log(x)))

.env$replaceNAs <- function(x) {
	nas <- is.na(x)
	if(sum(nas) > 0) {
		if (class(x) == "character") {
			x[nas] <- "<BLANK>"
		} else if (class(x) == "factor") {
			x <- as.character(x)
			x[nas] <- "<BLANK>"
			x <- as.factor(x)
		} else if (class(x) %in% c("numeric", "integer")) {
			x[nas] <- 0
		} else if (class(x) == "logical") {
			x[nas] <- FALSE
		}
	}
	x
}

.env$object.sizes <- function(plot=FALSE,value=TRUE) {
    os <- rev(sort(sapply(ls(envir=.GlobalEnv), function (o) object.size(get(o)))))
    if (plot) {
    	if (length(os) > 50) {
			dotchart(os[50:1], main="Memory usage by object", xlab="Bytes")
		} else {
			dotchart(rev(os), main="Memory usage by object", xlab="Bytes")
		}
	}
	if (value)
		return(os)
}

#### calculate a running function such as a running average
# x = data to calcuate function on
# order.by = data to first order x by (optional)
# FN = function to run
# extend = repeat first and last valid calculation to the ends,
#            otherwise the first and last window/2 values are NA
# exclude.center = do not include the center point itself in the 
#            analysis, just it's neighbors
# multicore = run the function in parallel using the 'multicore' library
# ... = additional arguments for function FN
.env$running <- function(x, order.by=NULL, FN="mean", window=10, extend=TRUE, exclude.center=FALSE, multicore=FALSE, ...)
{
    if (window == 1) return(x)
    FN=match.fun(FN)
    len <- length(x)
    before <- window %/% 2
    after <- window - before - 1
    start <- before + 1
    end <- len - after
    y <- numeric(len)

    if (! missing(order.by)) {
        if (length(x) != length(order.by)) {    
            warning("length(x) != length(order.by)")
            return(x)
        }
        o <- order(order.by)
        x <- x[o]
        orig_order <- order( (1:len)[o] )
        reorder <- TRUE
    } else {
        reorder <- FALSE
    }

    before.to.after <- -before:after
    if (exclude.center) before.to.after <- before.to.after[ before.to.after != 0 ]

    if (multicore) {
        require('multicore')
        sets <- mclapply(start:end, function(i) x[i + before.to.after])
        y <- unlist(mclapply(sets, FN, ...)) 
    } else {
      sets <- lapply(start:end, function(i) x[i + before.to.after])
      y <- base:::simplify2array(lapply(sets, FN, ...), higher=TRUE) 
    }

    if (extend == TRUE) {
        y <- c(rep(y[1], before),y, rep(y[length(y)], after))
    } else {
        y <- c(rep(NA, before),y, rep(NA, after))
    }

    if (reorder) y <- y[orig_order]

    return(y)
}


##### buffered writing of a matrix to a file
# specify ":cb" as file to write to clipboard
.env$write.matrix <- function (x, file = "", sep = "\t", blocksize=4000, use.rownames=FALSE, use.colnames=TRUE, strip.whitespace=TRUE, header=NULL)  {

	if ("tbl" %in% class(x) || "data.table" %in% class(x) ) x <- as.data.frame(x)
	if (identical(rownames(x), as.character(seq_len(nrow(x))))) use.rownames=FALSE

	cn <- colnames(x)
	rn <- rownames(x)
   	x <- as.matrix(x)
   	
   	if (is.null(cn)) use.colnames <- FALSE
   	if (is.null(rn)) use.rowname <- FALSE

	is.pipe <- FALSE
	if (file == ":cb") {
		is.pipe <- TRUE
		file <- pipe('pbcopy', 'w')
	}

    # strip white space
	if ((is.logical(strip.whitespace) && strip.whitespace == TRUE ) || (is.numeric(strip.whitespace) && length(strip.whitespace) > 0)) {
		if (is.numeric(strip.whitespace)) {
			cols.to.strip <- strip.whitespace
		} else {
			cols.to.strip <- 1:ncol(x)
		}
		for (the.col in cols.to.strip) {
			x[,the.col] <- sub("^ +","", x[,the.col])
		}
	}
    
	# add header
	if (! is.null(header)) {
		cat(header, file=file, sep="\n")
		append <- TRUE
	} else {
		append <- FALSE
	}

    p <- ncol(x)
    seps <- c(rep(sep, p - 1), "\n")
    if (use.rownames) {
    	if (! is.null(cn) && use.colnames) cn <- c("ROW", cn)
    	seps <- c(sep, seps)
    }


    if (!missing(blocksize) && blocksize > 0) {
    	if (use.colnames) cat(cn, file = file, sep = seps, append = append)
        nlines <- 0
        nr <- nrow(x)
        while (nlines < nr) {
            nb <- min(blocksize, nr - nlines)
            block <- x[nlines + (1:nb), ]
			if (use.rownames) block <- cbind(rn[nlines+(1:nb)],block)
            cat(t(block), file = file, append = TRUE, sep = seps)
            nlines <- nlines + nb
        }
    } else if (use.rownames) {
    	if (use.colnames) {
	    	cat(c(cn, t(cbind(rn,x))), file = file, sep = seps, append = append)
	    } else {
	    	cat(t(cbind(rn,x)), file = file, sep = seps, append = append)
	    }
    } else {
    	if (use.colnames) {
	    	cat(c(cn, t(x)), file = file, sep = seps, append = append)
	    } else {
	    	cat(t(x), file = file, sep = seps, append = append)
		}
    }
    
    if (is.pipe) close(file)
}
.env$write.cb <- function(x, ...) { .env$write.matrix(x, file=":cb", use.rownames=TRUE, use.colnames=TRUE, ...) }

##### identify runs of values which meet a certain cutoff
.env$runs <- function(data, cutoff, frac, any=F, min) {
    if (length(data) == 0) return(NA)
    if (missing(cutoff)) {
        warning("cutoff must be defined")
        return(NA)
    }
    count = length(data[!is.na(data) && data >= cutoff])
    if (count == 0) return(0)
    if (any == TRUE) return(1)
    if (missing(frac) && missing(min)) {
        warning("frac or min must be defined")
        return(NA)
    }
    if (! missing(min)) {
        if (count >= min) return(T)
        else return(F)
    } else {
        if (count / length(data) >= frac) return(T)
        else return(F)
    }
}

##### calculate a running wilcox p.value (based on the function running)
.env$running_wilcox <- function(x, x2, window=10, skip=3, extend=TRUE, ...)
{
    if (window == 1) return(x)
    len <- length(x)
    before <- window %/% 2
    after <- window - before - 1
    start <- before + 1
    end <- len - after
    y <- numeric(len)

    positions <- seq(start,end,skip)
    if (positions[length(positions)] != end) positions <- c(positions, end)
    last <- start

    for (i in positions) {
             z <- wilcox.test(x[(i-before):(i+after)], x2, ...)
             y[last:i] <- z$p.value
             last <- i + 1
    }
    
    if (extend == TRUE) {
        y[1:(start-1)] <- rep(y[start], before)
        y[(end+1):len] <- rep(y[end], after)
    } else {
        y[1:(start-1)] <- rep(NA, before)
        y[(end+1):len] <- rep(NA, after)
    }
    return(y)
}


#### identify local maxima in a moving window
.env$local_peak_simple <- function(x, min=4, window=10, single_peak=TRUE)
{
    if (window == 1) return(x)
    len <- length(x)
    before <- window %/% 2
    after <- window - before - 1
    start <- before + 1
    end <- len - after
    y <- x
    dummy <- 1:window
 
    i <- start
    while (i <= end) {
        ibefore <- i-before; iafter <- i+after
        # skip this range if center is NA
#       if (is.na(y[i])) next
        # save the range
        this_window <- y[ibefore:iafter]
        # if the whole range is NA, skip
        if (length(dummy[!is.na(this_window)]) == 0) {
            # increment i 
            i <- i + 1
            # go on
            next
        }
    
        # get the max for the window    
        window_max <- max(this_window,na.rm=T)
        
        # if it is below the min
        if (window_max < min) {
            # null out entire window
            this_window[] <- NA
            # move over to the first non null (but not too far)
            new_i <- iafter + 1
            if (new_i > end && i < end)
                i <- end
            else if (i < end)
                i <- new_i
    
        # otherwise, only keep the max for the window
        } else {
            # T/F for positions of maximums
            maxs <- !is.na(this_window) && this_window == window_max
            # position of first max
            first_max <- (dummy[maxs])[[1]]
            # null out all non maximums
            this_window[! maxs] <- NA
            # if we only want a single peak for the window
            if (single_peak == T && first_max < window)
                # null out everything after the first max
                this_window[(first_max+1):window] <- NA
            # move i over so that left end of next window
            # is at on first_max (skipping over any nulls)
            new_i <- i + first_max - 1
            if (first_max > i && new_i > end && i < end)
                i <- end
            else if (first_max > i)
                i <- new_i
            else 
                i <- i + 1
        }
        # save the window
        y[ibefore:iafter] <- this_window
    }
    
    return( ! is.na(y))
}

##### more sophisticated local maxima detection
.env$peaksAndValleys <- function(d, delta=diff(range(d))/100) {
	################## initial feature finding
	features <- list()
	cumdiff <- 0
	for (i in 2:length(d)) {
		cumdiff <- d[i] - d[i-1] + cumdiff
		if (cumdiff < -delta) {
			# check if we keep going down
			if (i < length(d) & d[i] > d[i+1]) next
			# otherwise it's a valley
			features[[ length(features) + 1]] <- data.frame(i=i, type="V")
			cumdiff <- 0
		} else if (cumdiff > delta) {
			# check if we keep going up
			if (i < length(d) & d[i] < d[i+1]) next
			# otherwise it's a peak
			features[[ length(features) + 1]] <- data.frame(i=i, type="P")
			cumdiff <- 0
		}
	}
	features <- do.call("rbind",features)


	########## add peaks and valleys at the ends
	#### check before first feature
	if (features$i[1] > 1) {
		before.first <- d[ 1:( features$i[1] - 1) ]
		# any higher points before the first feature?
		if (any(before.first > d[ features$i[1] ])) {
			# skip any very tiny peak before a first valley
			if (! (features$type[1] == "V" & max(before.first) < d[ features$i[1] ] + delta) )
				features <- rbind( data.frame(i=which.max(before.first), type="P"), features)
		}
		# any lower points before the first feature?
		if (any(before.first < d[ features$i[1] ])) {
			# skip any very tiny valley before a first peak
			if (! (features$type[1] == "P" & min(before.first) > d[ features$i[1] ] - delta) )
				features <- rbind( data.frame(i=which.min(before.first), type="V"), features)
		}
	}

	#### check after last feature
	lft <- nrow(features)
	ld <- length(d)
	if (features$i[lft] < ld) {
		after.end <- d[ (features$i[lft]+1):ld ]
		# any higher points after the last feature?
		if (any(after.end > d[ features$i[lft] ])) {
			# skip any very tiny peak after a last valley
			if (! (features$type[lft] == "V" & max(after.end) < d[ features$i[lft] ] + delta) ) 
				features <- rbind(features, data.frame(i=which.max(after.end) + features$i[lft], type="P"))
		}
		# any lower points after the last feature?
		lft <- nrow(features)
		if (any(after.end < d[ features$i[lft] ])) {
			# skip any very tiny valley after a last peak
			if (! (features$type[lft] == "P" & min(after.end) > d[ features$i[lft] ] - delta) )
				features <- rbind( features, data.frame(i=which.min(after.end) + features$i[lft], type="V"))
		}
		features <- features[ order(features$i), ]
	}

	################ collapse runs of peaks or valleys
	# simplify runs of peaks or valleys
	runs <- rle(features$type)$lengths
	for (i in which(runs > 1)) {
		from <- sum(runs[1:(i-1) ])+1
		to <- sum(runs[1:i])
		if (features$type[from] == "P") {
			# get rid of everything except for the highest peak
			toRemove <- (from:to)[-which.max( d[ features$i[from:to] ] ) ]
		} else if (features$type[from] == "V") {
			# get rid of everything except for the lowest valley
			toRemove <- (from:to)[-which.min( d[ features$i[from:to] ] ) ]
		}
		features$type[toRemove] <- "REMOVE"
	}
	features <- features[ features$type != "REMOVE", ]

	return(features)
}





## use rle for these:
##### return positions of first items in runs
# .env$just_first <- function(x) {
#     if (length(x) < 1) return(x)
#     prev <- x[1] - 1
#     x1 <- vector(mode="logical", length=length(x))
#     for (i in 1:length(x)) {
#         if (x[i] == prev) {
#             x1[i] <- F
#         } else {
#             x1[i] <- T
#             prev <- x[i]
#         }
#     }
#     return(x1)
# }
# 
# 
# ##### return positions of items in runs whenever there is a transition
# .env$just_transitions <- function(x) {
#     if (length(x) < 1) return(x)
#     x1 <- rep(F,length(x))
#     maxi <- length(x) - 1
#     i <- 2
#     while (1) {
#         if (i > maxi) break 
#         if (x[i] != x[i+1]) {
#             x1[c(i,i+1)] <- T
#             i <- i + 2
#         } else if (x[i] != x[i-1]) {
#             x1[i] <- T
#             i <- i + 1
#         } else {
#             i <- i + 1
#         }
#     }
#     x1[1] <- T
#     x1[maxi+1] <- T
#     return(x1)
# }

 
##### for use with snow parallel computing library
##### cl is a snow cluster object
# .env$parallel_running <- function(x, cl, FN="mean", window=10, extend=TRUE, ...)
# {
#     if (window == 1) return(x)
#     FN=match.fun(FN)
#     len <- length(x)
#     num_window <- len-window+1
# 
#     tmp_data <- x[1:num_window]
#     for (i in 2:window) tmp_data <- c(tmp_data, x[i:(num_window+i-1)])
#     windowed_array <- matrix(tmp_data,nrow=num_window,ncol=window,byrow=F)
#     rm(tmp_data)
# 
#     #### do all the work in one line:
#     y <- parRapply(cl,windowed_array, FN)
# 
#     before <- window %/% 2
#     after <- window - before - 1
#     if (extend) {
#         y <- c(rep(y[1], before),y,rep(y[length(y)],after))
#     } else {
#         y <- c(rep(NA, before),y,rep(NA,after))
#     }
#     return(y)
# }
# 
# .env$db.query <- function(group,query) {
#     library(RMySQL)
#     con <- dbConnect(dbDriver("MySQL"), group=group)
#     rs <- dbSendQuery(con,query)
#     d <- fetch(rs, -1)
#     dbDisconnect(con)
#     return(d)
# }
#   
# .env$TMAX.query <- function(q) {
#   return(.env$db.query("tmax",q))
# }
# 
# .env$CGH.query <- function(q) {
#   return(.env$db.query("CGH",q))
# }
# 
# .env$NCBI.query <- function(q) {
#   return(.env$db.query("ncbi_gene",q))
# }
# 
# .env$UCSC.query <- function(q) {
#   return(.env$db.query("ucsc",q))
# }
# 
# 
# .env$oracle_query <- function(q) {
# 	Sys.putenv(ORACLE_HOME="/usr/local/oracle/orahome")
# 	Sys.putenv(DYLD_LIBRARY_PATH="/usr/local/oracle/orahome/lib:/usr/local/oracle/orahome/rdbms/lib")
# 	Sys.putenv(CLASSPATH="/usr/local/oracle/orahome/jdbc/lib/classes12.zip")
#     library(ROracle)
#     con <- dbConnect(Oracle(), user="gcutler", password="gcutler", db="or028p")
#     rs <- dbSendQuery(con,q)
#     d <- fetch(rs, -1)
#     dbDisconnect(con)
#     return(d)
# }



.env$roc <- function(tp, metrics, do.plot=TRUE, round=NULL, max.thresholds=NULL) {

	roc.points <- lapply(data.frame(metrics), function(x) {
		are.na <- is.na(tp) || is.na(x)
		if (is.null(round)) {
			thresholds <- unique(sort(x[!are.na]))
		} else {
			thresholds <- unique(round(sort(x[!are.na]),round))
		}
		if (is.null(max.thresholds) || max.thresholds >= length(thresholds)) {
			thresholds <- c( thresholds[1] - 1, thresholds, thresholds[length(thresholds)]+1 )
		} else {
			thresholds <- c( thresholds[1] - 1, thresholds[ round(seq(1,length(thresholds),length=max.thresholds)) ],   thresholds[length(thresholds)]+1 )
		}
		mt <- matrix(thresholds)
		score.t <- apply(mt, 1, gt, values=x[ !are.na && tp == TRUE ])
		score.f <- apply(mt, 1, gt, values=x[ !are.na && tp == FALSE ])
		scores <- data.frame(tp=score.t, fp=score.f, threshold=thresholds)
		scores <- scores[ order(scores$tp), ]
		
		# trapezoid function
		idx = 2:length(score.f)
		auc <- (diff(scores$fp) %*% (scores$tp[idx] + scores$tp[idx - 1]))/2
		list(scores=scores , AUC=auc)
		
	})
	
	if (do.plot) {
		plot(0,0,xlim=c(0,1),ylim=c(0,1),xlab="false positive",ylab="true positive",main="ROC Curves",type="n")
		abline(a=0,b=1,lty=2,lwd=3,col="grey")
		colors <- rainbow(length(roc.points)+2)[1:length(roc.points)]
		for (i in 1:length(roc.points)) points(roc.points[[i]]$scores$tp ~ roc.points[[i]]$scores$fp, col=colors[i], type="o")
		legend("bottomright", paste(colnames(metrics), sapply(roc.points,function(x) round(x$AUC,2)), sep="="), col=colors, lwd=3)
	}

	roc.points
}

.env$AUC.trapezoid <- function(y, x=NULL) {
	if (is.null(x)) {
		if (length(y) >= 3) {
			1/2 * (sum(y) + sum(y[2:(length(y)-1)]))
		} else if (length(y) == 2) {
			(y[1] + y[2]) / 2
		} else {
			NA
		}
	} else if (length(x) != length(y)) {
		stop("length(y) and length(x) differ")
	} else if (length(y) < 2) {
		NA
	} else {
		sum( sapply(seq_along(y)[-1], function(i) (y[i-1]+y[i]) * (x[i]-x[i-1]) / 2 ) )
	}
}



.env$lt.diff <- function(threshold, values.ref=values.ref, values.exp=values.exp) {
    a <- sum(values.ref < threshold)
    b <- sum(values.exp < threshold)
    return(b/(a+b))
}


#### chromosome definitions
.env$chr.human <- data.frame(
	chrom=1:25,
    end=c(249250621,243199373,202541808,191412218,180915260,171115067,159138663,146364022,141213431,141948456,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,68956809,63025520,48129895,51304566,155270560,59373566,16571),
	global=c(0,249250621,492449994,694991802,886404020,1067319280,1238434347,1397573010,1543937032,1685150463,1827098919,1962105435,2095957330,2211127208,2318476748,2421008140,2511362893,2592558103,2670635351,2739592160,2802617680,2850747575,2902052141,3057322701,3116696267),
	colors=c("#996600","#666600","#99991E","#CC0000","#FF0000","#FF00CC","#FFCCCC","#FF9900","#FFCC00","#FFFF00","#CCFF00","#00FF00","#358000","#0000CC","#6699FF","#99CCFF","#00FFFF","#CCFFFF","#9900CC","#CC33FF","#CC99FF","#666666","#999999","#CCCCCC","#aaaaaa"),
    name=as.character(c(1:22,"X","Y","M"))
)
.env$chr.human$name <- as.character(.env$chr.human$name)

.env$chr.mouse <- data.frame(
       chrom=1:22,
       end=c(197062826,182491847,159871126,155175021,152002774,149521521,145133449,132071252,123986862,129936888,122076494,120458182,120614378,123956386,103482684,98373251,95170126,90727465,61319255,165555972,2724008,11542),
       global=c(0,197062826,379554673,539425799,694600820,846603594,996125115,1141258564,1273329816,1397316678,1527253566,1649330060,1769788242,1890402620,2014359006,2117841690,2216214941,2311385067,2402112532,2463431787,2628987759,2631711767),
       colors=c("#996600","#666600","#99991E","#CC0000","#FF0000","#FF00CC","#FFCCCC","#FF9900","#FFCC00","#FFFF00","#CCFF00","#00FF00","#358000","#0000CC","#6699FF","#99CCFF","#00FFFF","#CCFFFF","#9900CC","#CC33FF","#CC99FF","#CCCC99"),
       name=c(1:19,"X","Y","M")
)
.env$chr.mouse$name <- as.character(.env$chr.mouse$name)


.env$genome.plot <- function(chroms, yrange, xrange, chr.range, data, xlab, ylab, type, pch, label.chrs=T, zero.line=T, color="bychr", ...) {
	if (missing(chroms))
		chroms <- chr.human
	chrom.num <- nrow(chroms)
	plot.chroms <- 1:max(chroms$chrom)
	if (missing(xrange)) {
		if (! missing(chr.range)) {
		   xrange <- c(chroms$global[min(chr.range)], chroms$global[max(chr.range)+1])
		   plot.chroms <- min(chr.range):max(chr.range)
		} else if (! missing(data)) {
		   xrange <- range(data$x, na.rm=T)
		   min.chrom <- max(which(chroms$global < min(data$x, na.rm=T)))
		   if (length(min.chrom) == 0 || min.chrom < 1) min.chrom <- 1
		   max.chrom <- which(chroms$global > max(data$x, na.rm=T))
		   if (length(max.chrom) == 0) {
		   		max.chrom <- chroms$chrom[ nrow(chroms) ]
		   	} else if ( is.infinite(max.chrom) || max.chrom < 1 || max.chrom > chroms$chrom[ nrow(chroms) ]) {
		   		max.chrom <- chroms$chrom[ nrow(chroms) ]
		   	} else {
		   		max.chrom <- max.chrom[1]
		   	}
		   plot.chroms <- min.chrom:max.chrom
		} else {
		   xrange <- c(1,chroms$end[chrom.num]+chroms$global[chrom.num])
		}
	}

	if (missing(yrange) && ! missing(data)) {
		yrange <- range(data$y, na.rm=T)
	}
	
	if (missing(ylab))
		ylab <- "log2 Fold Change"
	if (missing(xlab))
		xlab <- "Global Genomic Position"
	plot(0, xlim=xrange,ylim=yrange,type="n", xlab=xlab, ylab=ylab, ...)
	
	for (i in plot.chroms) {
		abline(v=chroms$global[i], col="#996633")
		if (!is.null(label.chrs)) {
			if (class(label.chrs) == "logical" && label.chrs == TRUE) {
				text((chroms$global[i] + chroms$end[i]/2), yrange[2] * 1.025, chroms$name[i], col="darkgreen")
			} else if (class(label.chrs) == "numeric") {
				text((chroms$global[i] + chroms$end[i]/2), label.chrs, chroms$name[i], col="darkgreen")
			}
		}
	}
	abline(v=(chroms$end[max(plot.chroms)]+chroms$global[max(plot.chroms)]), col="#996633")
	if (zero.line)
		abline(h=0, col="grey", lwd=4)
	if (class(color) == "character") {
		if (color == "rainbow")
			colors <- rep(rainbow(nrow(chroms)+1))
		else if (color == "bychr") 
			colors <- chroms$colors
		chr.homes <- cut(data$x, c(chroms$global, chroms$global[nrow(chroms)]+chroms$end[nrow(chroms)]+ 1), labels=F)
		colors <- colors[ chr.homes ]
	} else if (class(color) == "numeric") {
		colors <- color
	} else if (is.null(color)) {
		colors <- 1
	}
		
		
	if (! missing(data)) {
		if (missing(type))
			type <- "p"
		if (missing(pch))
			pch <- 1
		points(data$y ~ data$x, type=type, pch=pch, col=colors)
	}
}

.env$genome_plot <- .env$genome.plot

.env$chi.square <- function(a,b) {
	contingency.table <- table(a,b)
	row.tot <- margin.table(contingency.table,1)
	col.tot <- margin.table(contingency.table,2)
	tot.tot <- sum(row.tot)
	
	ncol.cont.table <- ncol(contingency.table)
	nrow.cont.table <- nrow(contingency.table)
	
	row.tot.m <- matrix(rep(row.tot, ncol.cont.table), ncol=ncol.cont.table)
	col.tot.m <- matrix(rep(col.tot, each=nrow.cont.table), nrow=nrow.cont.table)
	
	expected <- col.tot.m / tot.tot * row.tot.m 
	
	contingency.squared.diff <- (contingency.table - expected)^2/expected
	chi.square <- sum(contingency.squared.diff)
	degree.freedom <- (ncol.cont.table-1) * (nrow.cont.table-1)
	p.val <- pchisq(chi.square,degree.freedom,lower=F)
	
	n <- length(a)
	k <- min(ncol.cont.table,nrow.cont.table)
	cramer <- sqrt(chi.square/ ( n * (k-1)))
	T.T.proportion.over.expected <- contingency.table[2,2] / expected[2,2]
	return( list(chi.square=chi.square, degrees.of.freedom=degree.freedom, p.value=p.val, cramers.phi=cramer, T.T.proportion.over.expected=T.T.proportion.over.expected))
}

.env$chi.square.quick <- function(a,b) {
	contingency.table <- table(a,b)
	row.tot <- margin.table(contingency.table,1)
	if (length(row.tot) < 2) return(NA)
	col.tot <- margin.table(contingency.table,2)
	if (length(col.tot) < 2) return(NA)
	tot.tot <- sum(row.tot)
	
	row.tot.m <- matrix(rep(row.tot, 2), ncol=2)
	col.tot.m <- matrix(rep(col.tot, each=2), nrow=2)
	
	expected <- col.tot.m / tot.tot * row.tot.m 
	
	contingency.squared.diff <- (contingency.table - expected)^2/expected
	return(pchisq(sum(contingency.squared.diff),1,lower=F))
}




.env$hypergeometric.analysis <- function(a, b) {
	# q = number of positives selected
	p.q <- sum(a && b) - 1
	# m = number of total positives
	p.m <- sum(a)
	# n = number of total negatives
	p.n <- length(a) - p.m
	# k = number of sampled items
	p.k <- sum(b)
	pval <- phyper(p.q, p.m, p.n, p.k, lower.tail=F)
	
	# calculate observed / expected	
	contingency.table <- table(a,b)
	row.tot <- margin.table(contingency.table,1)
	col.tot <- margin.table(contingency.table,2)
	tot.tot <- sum(row.tot)
	ncol.cont.table <- ncol(contingency.table)
	nrow.cont.table <- nrow(contingency.table)
	row.tot.m <- matrix(rep(row.tot, ncol.cont.table), ncol=ncol.cont.table)
	col.tot.m <- matrix(rep(col.tot, each=nrow.cont.table), nrow=nrow.cont.table)
	expected <- col.tot.m / tot.tot * row.tot.m 
	fold <- contingency.table[2,2] / expected[2,2]

	return( list(pval=pval, fold=fold) )
}

.env$count.na <- function(x) { sum(is.na(x)) }
.env$count.not.na <- function(x) { sum(! is.na(x)) }
.env$fraction.true <- function(x) { sum(x)/length(x) }
.env$original.order <- function(o) { order( (1:length(o))[o] ) }

.env$log.ratios <- function(x) log2(x)
.env$unlog.ratios <- function(x) x <- 2^x

.env$log.fold.change <- function(x) { .env$log.ratios(.env$fold.change.to.ratio(x)) }
.env$unlog.fold.change <- function(x) { .env$ratio.to.fold.change(.env$unlog.ratios(x)) }

.env$ratio.to.fold.change <- function(x) {
	x[ x == 0 ] <- NA
	na.x <- is.na(x)
	less.than.1 <- ! na.x && x < 1
	x[ less.than.1 ] <- -1 / x[ less.than.1 ]
	return(x)
}

.env$fold.change.to.ratio <- function(x) {
	x[ x == 0 ] <- NA
	na.x <- is.na(x)
	less.than.0 <- ! na.x && x < 0
	x[ less.than.0 ] <- -1 / x[ less.than.0 ]
	return(x)
}


.env$mouse.chromosome.to.int <- function(x) {
	x <- as.character(x)
	x[ x == "X" ] <- "20"
	x[ x == "Y" ] <- "21"
	x[ x == "M" ] <- "22"
	x[ x == "MT" ] <- "22"
	x <- as.integer(x)
	return(x)
}

.env$human.chromosome.to.int <- function(x) {
	x <- as.character(x)
	x <- sub("^chr","",x)
	x[ x == "X" ] <- "23"
	x[ x == "Y" ] <- "24"
	x[ x == "M" ] <- "25"
	x[ x == "MT" ] <- "25"
	x <- as.integer(x)
	return(x)
}

.env$local.to.global <- function(pos,chr,species="human") {
	if (class(pos) == "data.frame" && "POS" %in% colnames(pos) && "CHROM" %in% colnames(pos)) {
		chr <- pos$CHROM
		pos <- pos$POS
	}
	if (species == "human") {
		return( pos + .env$chr.human$global[ .env$human.chromosome.to.int(chr) ])
	} else if (species == "mouse") {
		return( pos + .env$chr.mouse$global[ .env$mouse.chromosome.to.int(chr) ])
	}
}

#### calculate a weighted running mean
.env$weighted.running.mean <- function(x, weights, exclude.center=TRUE, mask=NULL)
{
   	window <- length(weights)
    if (window == 1) return(x)
    if (window %% 2 == 0) {
    	center <- window %/% 2
    	weights <- weights[c(0:center,mean(weights[c(center,center+1)],na.rm=T),(center+1):window)]
    	window <- length(weights)
    }

   	window.center <-  window %/% 2 + 1 
    if (exclude.center) {
    	weights[ window.center] <- 0
    }
 
    weights <- weights / sum(weights)
    
   	do.mask <- F
    if (!is.null(mask)) {
	   	mask.weights <- rep(1,length(x))
    	mask.weights[ mask ] <- 0
    	do.mask <- T
	}

	len <- length(x)
    before <- window %/% 2
    after <- window - before - 1
    start <- before + 1
    end <- len - after
    y <- numeric(len)
    for (i in start:end) {
    	this.window <- (i-before):(i+after)
    	if (do.mask) {
    		this.mask <- mask.weights[this.window]
    		if (!exclude.center) this.mask[window.center] <- 1
	    	these.weights <- weights * this.mask
	    	these.weights <- these.weights/sum(these.weights)
			y[i] <- sum(x[this.window] * these.weights, na.rm=T)
		} else {
			y[i] <- sum(x[this.window] * weights, na.rm=T)
		}
    }
    
    #### lm for edges
    left.window <- 1:start
    right.window <- end:len
    if (do.mask) {
    	left.window <- left.window[ mask.weights[ left.window ] > 0 ]
    	right.window <- right.window[ mask.weights[ right.window ] > 0 ]
    }
    ## left edge
    lm.left <- lm( x[left.window] ~ I(1:length(left.window)), na.action=na.omit)
    y[left.window] <- predict(lm.left)[1:length(left.window)]
    ## right edge
    lm.right <- lm( x[right.window] ~ I(1:length(right.window)), na.action=na.omit)
    y[right.window] <- predict(lm.right)[1:length(right.window)]
    
    return(y)
}

##### adjust b towards a using a weighted neighborhood and differential mixing
.env$weighted.adjust <- function(a, b, min.weight=0)
{
	square.diff <- (b-a)^2
	mean.diff <- mean(square.diff, na.rm=T)
	sd.diff <- sd(square.diff,na.rm=T)
	z.diff <- square.diff/sd.diff
	z.diff <- log(z.diff - min(z.diff) + 1)
	adjust.w <- z.diff / max(z.diff) * (1 - min.weight)
	recip.adjust <- 1-adjust.w 

	len <- length(a)
    y <- numeric(len)
    for (i in 1:len) {
            y[i] <- (b[i] * adjust.w[i]) + (a[i] * recip.adjust[i])
    }
 	
    return(y)
}


.env$adjacent.centers <- function (x) {
	x1 <- x[-1]
	x2 <- x[ -length(x)]
	return( (x1+x2)/2 )
}

.env$truncate.abs <- function(x, a) {
	a <- abs(a)
	x[ x > a ] <- a
	x[ x < -a ] <- -a
	return(x)
}

.env$pvals.lmer<-function(model, conf.level=.95, conf=TRUE) {
   co<-summary(model)@coefs
   se<-co[,2]
   t<-abs(co[,3])
   df<-nrow(model@frame)-nrow(co)
   p.val<-(1-pt(t, df))*2
   ci.l<-co[,1] - qt(conf.level+(1-conf.level)/2, df)*se
   ci.u<-co[,1] + qt(conf.level+(1-conf.level)/2, df)*se
   ci<-cbind(ci.l, ci.u)
   colnames(ci)<-c(paste((1-conf.level)/2*100, "%", sep=""),
                   paste((conf.level+(1-conf.level)/2)*100, "%",sep=""))
   if (conf) {
     return(cbind(co, df, p.val, ci))
   }
   else {
     return(cbind(co, df, p.val))
   }
}


##### color functions
.env$saturate <- function(colors, add=.1, multiply=1) {
	require('fBasics')
	colors2 <- rgb2hsv(t(.asRGB(colors)))
	colors2[2,] <- (colors2[2,] + add)*multiply
	colors2[colors2 > 1] <- 1
	return(hsv(colors2[1,], colors2[2,], colors2[3,]))
}

.env$darken <- function(colors, add=.1, multiply=1) {
	require('fBasics')
	colors2 <- rgb2hsv(t(.asRGB(colors)))
	colors2[3,] <- (colors2[3,] - add)/multiply
	colors2[colors2 > 1] <- 1
	return(hsv(colors2[1,], colors2[2,], colors2[3,]))
}

.env$make.transparent <- function(colors, alpha=.5) {
	colors2 <- col2rgb(colors,alpha=FALSE)
	alpha <- floor(alpha * 255)
	if (length(alpha) == ncol(colors2)) {
		colors2 <- rbind(colors2, alpha)
		return(apply(colors2, 2, function(x) rgb(x[1],x[2],x[3],x[4], maxColorValue=255)))
	} else {
		return(apply(colors2, 2, function(x) rgb(x[1],x[2],x[3],alpha, maxColorValue=255)))
	}
}

## normal kernel weights
.env$normal.weights <- function(weight.length) {

	# how many points in the normal distribution do we
	# need to get the final desired length?

	# this formula was calculated to give roughly the
	# right number
	pred.len <- round(4.57 * weight.length + 1.345,0)
	if (pred.len %% 2 == 0) pred.len <- pred.len + 1


	# calculate the weights
	weights <- dnorm( seq(-10,10,length=pred.len), mean=0, sd=1)
	weights <- weights / max(weights)
	weights <- weights[ weights > 0.1 ]

	# remove any extra points if there are too many
	while (length(weights) > weight.length) {
		weights <- weights[-1]
		if (length(weights) > weight.length) {
			weights <- weights[ -length(weights) ]
		}
	}

	# tack on any extra linearly interpolated points if the
	# final number of points is too few
	while (length(weights) < weight.length) {
		weights <- c( weights[1] - (weights[2]-weights[1])/2, weights)
		if (weights[1] < 0.01) weights[1] <-  0.01
		if (length(weights) < weight.length) {
			weights <- c( weights, weights[length(weights)] - (weights[length(weights)-1]-weights[length(weights)])/2)
			if (weights[length(weights)] < 0.01) weights[length(weights)] <-  0.01
		}
	}
	return(weights)
}

.env$running.weighted.mean.sd <- function(x, weights, exclude.center=F)
{
	window <- length(weights)
	weights <- weights / min(weights[weights > 0])
    len <- length(x)
	the.matrix <- matrix(0,ncol=window,nrow=(len-window+1))
	for (i in 1:window) the.matrix[,i] <- x[(window-i+1):(len-i+1)]

	if (exclude.center) {
		the.matrix <- the.matrix[, -(window %/% 2)]
		weights <- weights[ -(window %/% 2)]
	}

	means <- apply(the.matrix, 1, function(z) mean(rep(z,weights),na.rm=T))
	sds <- apply(the.matrix, 1, function(z) sd(rep(z,weights),na.rm=T))

	before <- window %/% 2
	after <- window - before - 1

	means <- c(rep(means[1], before),  means, rep(means[length(means)],after))
	sds <- c(rep(sds[1], before),  sds, rep(sds[length(sds)],after))
    return(list(mean=means,sd=sds))
}

.env$running.mean.sd <- function(x, window=10, exclude.center=F)
{
    len <- length(x)
	the.matrix <- matrix(0,ncol=window,nrow=(len-window+1))
	for (i in 1:window) the.matrix[,i] <- x[(window-i+1):(len-i+1)]

	if (exclude.center) {
		the.matrix <- the.matrix[, -(window %/% 2)]
	}

	means <- apply(the.matrix, 1, function(z) mean(z,na.rm=T))
	sds <- apply(the.matrix, 1, function(z) sd(z,na.rm=T))

	before <- window %/% 2
	after <- window - before - 1

	means <- c(rep(means[1], before),  means, rep(means[length(means)],after))
	sds <- c(rep(sds[1], before),  sds, rep(sds[length(sds)],after))
    return(list(mean=means,sd=sds))
}


.env$trimmed.running.weighted.mean.sd <- function(x, weights, trim)
{
	window <- length(weights) + 2*trim
	weights <- weights / min(weights[weights > 0])
    len <- length(x)
	the.matrix <- matrix(0,ncol=window,nrow=(len-window+1))
	for (i in 1:window)
		the.matrix[,i] <- x[(window-i+1):(len-i+1)]

	orders <- apply(the.matrix,1, order)
	to.na <- c(1:trim, (window-trim+1):window)
	orders <- t(apply(orders,2, function(z) {z[ z %in% to.na] <- NA; z}))
	the.matrix[ is.na(orders) ] <- NA

	means <- apply(the.matrix, 1, function(z) mean(rep(not.na(z),weights),na.rm=T))
	sds <- apply(the.matrix, 1, function(z) sd(rep(not.na(z),weights),na.rm=T))

	before <- window %/% 2
	after <- window - before - 1

	means <- c(rep(means[1], before),  means, rep(means[length(means)],after))
	sds <- c(rep(sds[1], before),  sds, rep(sds[length(sds)],after))
    return(list(mean=means,sd=sds))
}


.env$qq <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
    # Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
    if (class(pvector) == "data.frame") {
    	pvector <- as.matrix(pvector)
    	dim(pvector) <- NULL
    } else if (class(pvector) == "matrix") {
    	dim(pvector) <- NULL
	}    
    o <- -log10(sort(pvector,decreasing=F))
    N <- length(o)
    e <- -log10( 1:N/N )

	## create the confidence intervals
	c95 <- rep(0,N)
	c05 <- rep(0,N)
	
	## the jth order statistic from a 
	## uniform(0,1) sample 
	## has a beta(j,n-j+1) distribution 
	## (Casella & Berger, 2002, 
	## 2nd edition, pg 230, Duxbury)
	
	for(i in 1:N){
		c95[i] <- qbeta(0.95,i,N-i+1)
		c05[i] <- qbeta(0.05,i,N-i+1)
	}
	c95 <- -log10(c95)
	c05 <- -log10(c05)
	conf.int.points <- c(c95, rev(c05))
	conf.int.points[ conf.int.points > max(e) ] <- max(e)

	plot(0, ylim=c(0,max(e)), xlim=c(0,max(e)), type="n", 
	xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),
	main=title)


    polygon(c(e,rev(e)), conf.int.points, border=NA, col="#cfcfcf")
    lines(e,e,col="red")
    points(e,o,pch=19,cex=0.25)
	


    #You'll need ggplot2 installed to do the rest
#    plot= qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + 
#	    stat_abline(intercept=0,slope=1, col="red", xlim=c(0,max(e)), ylim=c(0,max(o)))
#    plot=plot+opts(title=title)
#    plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
#    plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
#    if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
#    plot
}


# test for proportions
# what fraction is test.set of reference sets
.env$proptest <- function(test.set, reference.set, full.set.length, alternative="greater") {
	if (class(reference.set) == "list") {
		pt.result <- lapply(reference.set, function(x) prop.test(
			c( sum(x %in% test.set), length(test.set)),c( length(x), full.set.length),
			alternative=alternative
		))
		sapply(pt.result, function(x) list(p.value=x$p.value,  ratio=x$estimate[1]/x$estimate[2]))
	} else {
		pt.result <- prop.test(
			c( sum(reference.set %in% test.set), length(test.set)), c( length(reference.set), full.set.length),
			alternative=alternative
		)
		list(p.value=pt.result$p.value, ratio=pt.result$estimate[1]/pt.result$estimate[2])
	}
}
.env$quick.anova.lm.pval <- function (formula, which.coefficient=1) {
  #   mf <- match.call(lm, call("lm", yg[5,] ~ cltg, i)) # for debuging outside of function
      mf <- match.call(expand.dots = FALSE)
          mf[[1L]] <- as.name("model.frame")
          mf <- mf[-3]
          mf <- eval(mf, parent.frame())
          y <- model.response(mf, "numeric")
          x <- model.matrix(attr(mf, "terms"), mf, NULL)
          np <- dim(x)
          n <- np[1]
          pX <- np[2]
      # original with dqrls doesn't work anymore as of R 3.0
      #    z <- .Fortran("dqrls", qr = x, n = n, pX = pX, y = y, ny = 1L, tol = 1e-7, coefficients = numeric(pX), residuals = y, effects = y, rank = integer(1), pivot = 1L:pX, qraux = double(pX), work = double(2 * pX), PACKAGE = "base")
      z <- .Call(stats:::C_Cdqrls, x, y, tol=1e-7)

          dfr <- n - z$rank
          p1 <- 1L:z$rank
          asgn <- attr(x, "assign")[ z$pivot ][p1]
          df <- c(unlist(lapply(split(asgn, asgn), length)), dfr)

          ssr <- sum(z$residuals^2)
          ss <- c(unlist(lapply(split(z$effects[p1]^2, asgn), sum)), ssr)
          f <- (ss * dfr) / (df * ssr ) #### F
          P <- .Internal(pf(f, df, dfr, FALSE, FALSE))
          P[which.coefficient]
    }

.env$gt.to.simple <- function(gt) {
	simple <- rep("WT", length(gt))
	simple[ is.na(gt) ] <- "NC"
	simple[ c(grep("^0/[1-9]", gt), grep("^[1-9]+/0", gt)) ] <- "HET"
	simple[ grep("^[1-9]/", gt) ] <- "VAR"
	simple[ gt == "0/0" ] <- "WT"
	simple[ gt %in% c(".","./.","") ] <- "NC"
	simple
}



.env$expandFactorColumns <- function(the.data, data.name="") {
	require('caret')
	wanted.cols <- colnames(the.data)
	wanted.cols <- wanted.cols[
				wanted.cols %nin% c("i","CHROM","POS","ID","REF","ALT", grep("\\.NUCLEOTIDE\\.ALT$", wanted.cols, value=T), grep("\\.GT$", wanted.cols,value=T)) &
				sapply(the.data, class) %nin% c("logical","integer","numeric")
					]
	for (the.col in wanted.cols) {
		i <- which(colnames(the.data) == the.col)
		if (length(unique(the.data[,i])) > 100) {
			warning("Too many levels for ", the.col,". Skipping", call.=FALSE)
			next
		# hold on to gt_simple column for later checking
		} else if (length(grep("\\.GT_SIMPLE",the.col)) == 1) {
			gtsimple <- the.data[,i, drop=FALSE]
		}

		dv <- dummyVars( paste("~ ", the.col, sep=""), data=the.data, sep=".")
		if (! is.null(dv)) {
			fm <- predict(dv , newdata=the.data)
#			fm <- data.frame(lapply(fm, as.factor))
			cat("dummyVariables for",data.name,the.col,":",ncol(fm),"\n")
			if (i > 1 && i < ncol(the.data)) {
				the.data <- data.frame(the.data[,1:(i-1), drop=FALSE], fm, the.data[,(i+1):ncol(the.data), drop=FALSE])
			} else if (i == 1) {
				the.data <- data.frame(fm, the.data[,-1, drop=FALSE])
			} else {
				the.data <- data.frame(the.data[,-ncol(the.data), drop=FALSE], fm )
			}
			rm(fm)
		}
		rm(dv,i)
		if (exists("gtsimple")) {
			the.data <- data.frame(the.data,gtsimple)
			rm(gtsimple)
		}
	}
	return(the.data)
}

.env$getLowInfoColumns <- function(x, min.frac=0.5) {
	q <- length(unique(x))
	if (q == 1) {
		return(TRUE)
	} else if (q > 1000) {
		return(FALSE)
	} else {
		q <- sort(table(x), decreasing=TRUE)
		if (q[2] < min.frac * length(x)) return(TRUE)
	}
	return(FALSE)
}

.env$isGzipped <- function(file.name) { length(grep("\\.gz$", file.name)) == 1 }

.env$VCFHeaderLineCount <- function(vcf.file, return.header=TRUE) {
	vcf.file.is.gzipped <- isGzipped(vcf.file)
	if (vcf.file.is.gzipped) {
		connection <- gzfile(vcf.file)
		tmp.txt <- scan(file=connection, what='character', sep="\n", n=300)
		close(connection)
	} else {
		tmp.txt <- scan(file=vcf.file, what='character', sep="\n", n=300)
	}
	lines.to.skip <- grep("^#CHROM\tPOS", tmp.txt)
	if (length(lines.to.skip) != 1 && lines.to.skip > 1) {
		warning(vcf.file, " is malformed\n", call.=FALSE)
		q(save="no", status = 1)
	} else {
		lines.to.skip <- lines.to.skip - 1
		header <- tmp.txt[1:lines.to.skip]
		rm(tmp.txt)
		cat(lines.to.skip, "header lines in vcf\n")
	}
	if (return.header) {
		return(list(lines=lines.to.skip, header=header))
	} else {
		return(lines.to.skip)
	}
}

.env$readVCF <- function(vcf.file) {
	header.info <- VCFHeaderLineCount(vcf.file)
	gz <- isGzipped(vcf.file)

	if (gz) {
		vcf <- read.delim(file=gzfile(vcf.file), skip=header.info$lines)
	} else {
		vcf <- read.delim(file=vcf.file, skip=header.info$lines)
	}
	colnames(vcf)[1] <- "CHROM"

	return(list(vcf=vcf, header=header.info$header))
}


####### machine learning functions


# TPR = (ref.event & pred.event) / ref.event = fraction of all events are predicted
# FPR = (ref.no_event & pred.event) / ref.no_event  = fraction of all non-events of are predicted
# PPV = (ref.event & pred.event) / pred.event	= fraction of predicted events that are true
# NPV = (ref.no_event & pred.no_event) / pred.no_event  = fraction of predicted non-events that are true
# ACC = (ref.event & pred.event + ref.no_event & pred.no_event) / (pred.event + pred.no_event)  = fraction of all calls that are true (var or wt)

.env$accuracyMetrics <- function(test,ref) {
	test <- test!="WT" && test!="NC"
	ref <- ref!="WT" && ref!="NC"

	data.frame(
		TPR = sum(ref & test)/sum(ref),
		FPR = sum(!ref & test)/sum(!ref),
		PPV = sum(ref & test)/sum(test),
		NPV = sum(!ref & !test)/sum(!test),
		ACC = (sum(ref & test) + sum(!ref & !test))/length(ref)
		)
}


.env$waitForUser <- function(waitFor=NULL) {
	if (is.null(waitFor)) {
		cat("Press [enter] to continue\n")
		flush.console()
		readline()
	} else {
		cat(paste0("Enter [", waitFor, "] to continue\n"))
		flush.console()
		while (1) {
			user.line <- readline()
			if (length(grep(waitFor,user.line))>0) break
		}
	}
	return(TRUE)
}


.env$splitIntoSets <- function(x, numberPerSet=1) {
	if (is.null(x)) stop("No argument supplied")
	if (length(x) <= 1) stop("Argument is empty")
	split( x, ceiling( 1:(length(x)) / numberPerSet ))
	}

.env$truncRange <- function(x, low=NULL, high=NULL, filter=FALSE) { 
	if (filter) {
		if (!is.null(low)) x <- x[x>=low]
		if (!is.null(high)) x <- x[x<=high]
	} else {
		if (!is.null(low)) x[x<low] <- low
		if (!is.null(high)) x[x>high] <- high
	}
	return(x)
}

.env$setcomp <- function(x, y, return.values=FALSE, silent=FALSE) {
    just.x <- setdiff(x, y)
    just.y <- setdiff(y, x)
    both <- intersect(y, x)
    if (! silent) {
        if (length(just.x) == 0 & length(just.y) == 0) {
            writeLines("Both sets equal")
        } else if (length(just.x) == 0) {
            writeLines(paste0("All x in y\n", length(just.y)," unique in y"))
        } else if (length(just.y) == 0) {
            writeLines(paste0(length(just.x)," unique in x\nAll y in x" ))
        } else {
            writeLines(paste0(length(just.x)," unique in x\n",  length(just.y)," unique in y\n", length(both)," in both\n"))
        }
    }
    if (return.values) {
        return(list(justX=just.x,justY=just.y,both=intersect(x,y)))
    }
}

.env$rmse.train_test <- function(model, data, truthColumn, train, drawPlot, ...) {
	prediction <- predict(model, data)
	if (drawPlot) {
		plot(prediction[train] ~ data[[truthColumn]][train], pch=21, bg="#377EB8", col="black" , 
		ylim=range(prediction), xlim=range(data[[truthColumn]]), xlab="Truth", ylab="Prediction" )
		points(prediction[-train] ~ data[[truthColumn]][-train], pch=21, bg="#E41A1C", col="black" )
	}
	results <- c( RMSE.Train=sqrt( mean( (prediction[train] - data[[truthColumn]][train])^2, na.rm=TRUE)),
			RMSE.Test=sqrt( mean( (prediction[-train] - data[[truthColumn]][-train])^2, na.rm=TRUE)),
			RMSE.TrainRandom=sqrt( mean( (sample(prediction[train]) - data[[truthColumn]][train])^2, na.rm=TRUE)),
			RMSE.TestRandom=sqrt( mean( (sample(prediction[-train]) - data[[truthColumn]][-train])^2, na.rm=TRUE))
		)
	results["Test.vs.Train"] <- results['RMSE.Test'] / results['RMSE.Train']
	results
}

.env$accuracy <- function(truth, prediction, train=NULL, randomize=FALSE) {
	if (!is.null(train)) {
		vals <- data.table(T=truth, P=prediction, train=seq_along(truth) %in% train)
		if (randomize) {
			vals[train, P := sample(P) ]
			vals[!train, P := sample(P) ]
		}			
		accuracy.train <- nrow(vals[train & T==P]) / nrow(vals[train])
		accuracy.test <- nrow(vals[!train & T==P]) / nrow(vals[!train])
		return(c(train=accuracy.train, test=accuracy.test))
	} else {
		if (randomize) truth <- sample(truth)
		return(accuracy=sum(truth==prediction) / length(truth))
	}
}


.env$partitionData <- function(x, p=0.75, cuts=5) {
	DF <- data.frame(x=x,i=seq_along(x))
	DF <- DF[ order(DF$x), ]
	DF$cut <- floor(DF$i * (cuts/(nrow(DF)+1)))
	set <- sort(unlist(by(DF$i, DF$cut, function(vals) sample(vals, length(vals)*p))))
	names(set) <- NULL
	set
}


.env$varImpGenes <- function(fitObject, geneMap, n=50, ...) {
	imps <- varImp(fitObject, ...)$importance
	imps <- imps[ order(imps$Overall, decreasing=TRUE),, drop=FALSE ]
	m <- use.match(rownames(imps), geneMap$Name)
	imps$Gene <- ""
	imps$Gene[m$a.order] <- geneMap$GeneID[m$b.order]
	if (is.null(n)) return(imps)
	imps[1:n,]
}

.env$predictionCorrelation <- function(model, data, truthColumn, trainIndex=NULL) {
	if (any(class(data) == 'data.table')) data <- as.data.frame(data)
	if (is.numeric(trainIndex)) trainIndex <- 1:nrow(data) %in% trainIndex
	if (is.null(trainIndex)) {
		prediction <- data.table(A=data[[truthColumn]], B=predict(model, data))
		cor(prediction$A, prediction$B)
	} else {
		prediction <- data.table(A=data[[truthColumn]], B=predict(model, data), C=trainIndex)
		results <- c(with(prediction[prediction$C,], cor(A, B)), with(prediction[!prediction$C,], cor(A, B)))
		names(results) <- c("TRAIN","TEST")
		results
	}
}

.env$varImpDplyr <- function(x, ...) {
	vi <- varImp(x, ...)
	if (class(vi) == "varImp.train") vi <- vi$importance
	vi %>% 
	as.data.frame() %>% 
	rownames_to_column(var="Feature") %>% 
	as.tibble %>% 
	dplyr::select(Feature, Importance=Overall, everything()) %>% 
	dplyr::arrange(-Importance)
}

	

.env$varImpPlot <- function(fitObject, n=25, quiet=FALSE) {
	vi <- varImpDplyr(fitObject)
	vi <- vi %>% dplyr::filter(Importance > 0) %>% dfilter(1:n() <= n)
	if (!quiet) print(vi)
	vi <- vi %>% arrange(Importance) %>% mutate(x=seq_len(n()))
	vi %>%
		ggplot(aes(y=Importance, x=x)) + 
		geom_bar(stat='identity', fill="darkorange") +
		scale_x_continuous(breaks=vi$x, labels=vi$Feature) +
		coord_flip()  + xlab("Feature") + ylab("Importance") -> p
	return(p)
}



### adapt from ggbiplot
.env$ggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
    obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
    ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
    alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
    varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
	## GC change:
    point.size = 2, text.col="black", group.palette = NULL,
    ...) 
{
    library(ggplot2)
    library(scales)
    library(grid)
    stopifnot(length(choices) == 2)
    if (inherits(pcobj, "prcomp")) {
        nobs.factor <- sqrt(nrow(pcobj$x) - 1)
        d <- pcobj$sdev
        u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
        v <- pcobj$rotation
    }
    else if (inherits(pcobj, "princomp")) {
        nobs.factor <- sqrt(pcobj$n.obs)
        d <- pcobj$sdev
        u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
        v <- pcobj$loadings
    }
    else if (inherits(pcobj, "PCA")) {
        nobs.factor <- sqrt(nrow(pcobj$call$X))
        d <- unlist(sqrt(pcobj$eig)[1])
        u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
        v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
            1]), FUN = "/")
    }
    else if (inherits(pcobj, "lda")) {
        nobs.factor <- sqrt(pcobj$N)
        d <- pcobj$svd
        u <- predict(pcobj)$x/nobs.factor
        v <- pcobj$scaling
        d.total <- sum(d^2)
    }
    else {
        stop("Expected a object of class prcomp, princomp, PCA, or lda")
    }
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
        FUN = "*"))
    v <- sweep(v, 2, d^var.scale, FUN = "*")
    df.v <- as.data.frame(v[, choices])
    names(df.u) <- c("xvar", "yvar")
    names(df.v) <- names(df.u)
    if (pc.biplot) {
        df.u <- df.u * nobs.factor
    }
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
    v.scale <- rowSums(v^2)
    df.v <- r * df.v/sqrt(max(v.scale))
    if (obs.scale == 0) {
        u.axis.labs <- paste("standardized PC", choices, sep = "")
    }
    else {
        u.axis.labs <- paste("PC", choices, sep = "")
    }
    u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
        100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
    if (!is.null(labels)) {
        df.u$labels <- labels
    }
    if (!is.null(groups)) {
        df.u$groups <- groups
    }
    if (varname.abbrev) {
        df.v$varname <- abbreviate(rownames(v))
    }
    else {
        df.v$varname <- rownames(v)
    }
#	df.v$varname[ grepl("_$", df.v$varname) ] <- ""  ## GC CHANGE
    df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
    g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
        ylab(u.axis.labs[2]) + coord_equal()
    if (var.axes) {
        if (circle) {
            theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                length = 50))
            circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                sin(theta))
            g <- g + geom_path(data = circle, color = muted("white"), 
                size = 1/2, alpha = 1/3)
        }
        g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
            xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
            "picas")), color = text.col)   ## GC change
    }
    if (!is.null(df.u$labels)) {
        if (!is.null(df.u$groups)) {
            g <- g + geom_text(aes(label = labels, color = groups), 
                size = labels.size)
        }
        else {
            g <- g + geom_text(aes(label = labels), size = labels.size)
        }
    }
    else {
        if (!is.null(df.u$groups)) {
            g <- g + geom_point(aes(color = groups), alpha = alpha, size = point.size)   ## GC change
        }
        else {
            g <- g + geom_point(alpha = alpha, size = point.size)   ## GC change
        }
    }
    if (!is.null(df.u$groups) && ellipse) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
        circle <- cbind(cos(theta), sin(theta))
        ell <- ddply(df.u, "groups", function(x) {
            if (nrow(x) <= 2) {
                return(NULL)
            }
            sigma <- var(cbind(x$xvar, x$yvar))
            mu <- c(mean(x$xvar), mean(x$yvar))
            ed <- sqrt(qchisq(ellipse.prob, df = 2))
            data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                mu, FUN = "+"), groups = x$groups[1])
        })
        names(ell)[1:2] <- c("xvar", "yvar")
        g <- g + geom_path(data = ell, aes(color = groups, group = groups))
    }
    if (var.axes) {
        g <- g + geom_text(data = df.v, aes(label = varname, 
            x = xvar, y = yvar, angle = angle, hjust = hjust), 
            color = text.col, size = varname.size)   ## GC change
            
        if (length(g$layers) == 3) g$layers <- g$layers[ c(2,1,3) ]  ## GC change
    }
    return(g)
}

.env$enrichmentAUC <- function(geneArray, signalArray, geneSet, randomize=NULL, plot=FALSE, plotTitle=NULL, returnTable=FALSE) {
	if (is.null(geneArray) || is.null(signalArray) || is.null(geneSet)) stop("Error: please supply gene list, signal values, and gene set list")
	if (length(geneArray) != length(signalArray)) stop("Error: Gene list is not the same length as signal list")
	if (plot) require(ggplot2)
	n <- length(geneArray)
	nSelected <- length(geneSet)
	nNonSelected <- n - nSelected
	scoreSelected <-  1.0 / nSelected
	scoreNonSelected <-  1.0 / nNonSelected
	
	o <- order(-signalArray)
	geneArray <- geneArray[o]
	signalArray <- signalArray[o]
	zero.point <- which.min(abs(signalArray))

	inSet <- geneArray %in% geneSet
	geneScore <- numeric(n)
	geneScore[inSet] <- scoreSelected
	geneScore[ !inSet ] <-  -scoreNonSelected
	cumGeneScore <- cumsum(geneScore)
	orderedScores <- sort(geneScore,decreasing=TRUE)
	bestCumGeneScore <- cumsum(orderedScores)
	bestNegCumGeneScore <- cumsum(rev(orderedScores))

	## kernel is normal distrubtion with mean = 0 and stdev = zero.point / 4
	enrichment.kernel <- c( dnorm(seq(0,4,length=zero.point), mean=0, sd=1), rep(0.0, n - zero.point) )
	depletion.kernel <- c( rep(0.0, zero.point), rev(dnorm(seq(0,4,length=n-zero.point), mean=0, sd=1)))

	maxEnrich <- sum( bestCumGeneScore * enrichment.kernel )
	maxDeplete <- sum( bestNegCumGeneScore * depletion.kernel ) # should be negative

	enrichScore <- sum( cumGeneScore * enrichment.kernel ) / maxEnrich
	depleteScore <- sum( cumGeneScore * depletion.kernel ) / maxDeplete # becomes positive
	if (depleteScore > 0 & enrichScore < depleteScore) {
		combined <- -depleteScore
		direction <- "dep"
	} else {
		combined <- enrichScore
		direction <- "en"
	}

	if (returnTable || plot) plotDT <- data.table(gene=geneArray, signal=signalArray, genePosition=seq_len(n), enrichmentScore=cumGeneScore)
	if (plot) {
		p <- ggplot(plotDT, aes(x=genePosition, y=enrichmentScore)) + geom_hline( yintercept=0, linetype="dashed") +
			geom_line() + ggtitle( plotTitle )
		print(p)
	}

	if (!is.null(randomize)) {
		randScore <- lapply(seq_len(randomize), function(x) {
			cum <- cumsum(sample(geneScore, replace=FALSE))
			c(en=sum( cum * enrichment.kernel ), dep=sum( cum * depletion.kernel ))
		})
		randScore <- as.data.frame(do.call("rbind", randScore))
		randScore$en <- randScore$en/maxEnrich
		randScore$dep <- randScore$dep/maxDeplete
		randScore$both <- randScore[[direction]] * ifelse(direction=="en", 1.0, -1.0)
		out <- c(enrichment=enrichScore, depletion=depleteScore, combined=combined, 
			enrichP.val=sum(randScore$en >= enrichScore)/randomize, 
			depleteP.val=sum(randScore$dep >= depleteScore)/randomize,
			combinedP.val=ifelse(direction == "en", sum(randScore$both > combined)/randomize, sum(randScore$both < combined)/randomize)
			)
	} else {
		out <- c(enrichment=enrichScore, depletion=depleteScore, combined=combined, 
			enrichP.val=NA,
			depleteP.val=NA,
			combinedP.val=NA
			)
	}
	if (returnTable) {
		out <- as.list(out)
		out[['table']] <- plotDT
	}
	return(out)
}

.env$plotEnrichments <- function(x, labelLines=TRUE, returnPlot=TRUE, plotTitle=NULL) {
	require(ggplot2)
	if (labelLines) require(ggrepel)
	# x should be a matrix/data frame/data table with cumsum scores, one column for each gene set
	cnames <- colnames(x)
	nr <- nrow(x)
	nc <- ncol(x)
	x <- as.matrix(x)
	dim(x) <- NULL
	long <- data.frame(set=rep(cnames, each=nr), value=x, position=rep(seq_len(nr), nc))
	p <- ggplot(long, aes(x=position, y=value, color=set, fill=set)) +
		geom_hline( yintercept=0, linetype="dashed") +
		geom_line() + ggtitle( plotTitle )
	if (labelLines) {
		labelLocs <- long %>% mutate(absval = abs(value)) %>% group_by(set) %>% filter( absval==max(absval) ) %>% ungroup
		p <- p + scale_color_discrete(guide=FALSE) +
			geom_label_repel(data=labelLocs, aes(label=set), 
				color="white", segment.color="#000000", segment.size=1, 
				fontface = 'bold', size=2.5, box.padding = unit(0.25, "lines"), 
				point.padding = unit(0.5, "lines"))
	} else {
		p <- p + scale_color_discrete(name="Gene Set")
	}
	print(p)
	if (returnPlot) return(p)
}


.env$membership.matrix <- function(x) {
	m <- matrix(NA, ncol=length(x), nrow=length(x))
	colnames(m) <- names(x)
	rownames(m) <- names(x)
	for (i in seq_along(x)) {
		for (j in i:length(x)) {
			m[j,i] <- m[i,j] <- sum(x[[i]] %in% x[[j]])
		}
	}
	m
}

.env$xvert <- function() theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.95))
.env$noaxis <- function(the.axis='x') {
	require(ggplot2)
	eb <- element_blank()
	if (the.axis == 'y') {
		return(theme(axis.text.y=eb, axis.ticks.y=eb))
	} else if (the.axis == 'x') {
		return(theme(axis.text.x=eb, axis.ticks.x=eb))
	} else {
		return(theme(axis.text.y=eb, axis.ticks.y=eb, axis.text.x=eb, axis.ticks.x=eb))
	}
}
.env$nolegend <- function() theme(legend.position="none")
		
	
.env$ggmultiplot <- function(..., plotlist=NULL, file, cols=NA, layout=NULL) {
	require(grid)
	plots <- c(list(...), plotlist)
	numPlots = length(plots)
	if (is.null(layout)) {
		if (is.na(cols)) cols <- ceiling(sqrt(numPlots))
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
	}
	if (numPlots==1) {
		print(plots[[1]])
	} else {
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		for (i in 1:numPlots) {
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
		}
	}
}

.env$plateHeatmap <- function(data, columnName="col", rowName="row", valueName="value", logScale=FALSE,  return=FALSE) {
	require(data.table)
	require(ggplot2)
	x <- as.data.table(data)[, .SD, .SDcols=c(columnName, rowName, valueName)]
	setnames(x, c("col","row","value"))
	if (logScale) x[, value := log10(value)]
	p <- ggplot(x) + geom_raster(aes(x=col, y=row, fill=value))

	if (return) {
		return(p)
	} else {
		print(p)
	}
}

.env$dselect <- dplyr::select
.env$dfilter <- dplyr::filter

.env$fun_if <- function(x, FUN, condition) {
    FUN=match.fun(FUN)
	x[ condition ] <- FUN( x[condition] )
	x
}

.env$crossplot <- function(x, XVAR, YVAR, CLASS, xlab=NA, ylab=NA, return=FALSE) {
	require(ggplot2)
	require(ggrepel)
	require(data.table)
	
	x2 <- data.table(x)[, .SD, .SDcols=c(XVAR,YVAR, CLASS)]
	setnames(x2, c(XVAR,YVAR,CLASS), c("XVAR","YVAR","CLASS"))
	if (is.na(xlab)) xlab <-  XVAR
	if (is.na(ylab)) ylab <-  YVAR

	p <- ggplot(x2[, .(
		medX=0.0+median(XVAR), XVAR.Q1=0.0+quantile(XVAR, probs=0.25), XVAR.Q3=0.0+quantile(XVAR, probs=0.75),
		medY=0.0+median(YVAR), YVAR.Q1=0.0+quantile(YVAR, probs=0.25), YVAR.Q3=0.0+quantile(YVAR, probs=0.75)
		), by=CLASS], aes(x=medX, y=medY, color=CLASS)) +
		geom_segment(aes(x=XVAR.Q1, xend=XVAR.Q3, y=medY, yend=medY, color=CLASS), size=3, alpha=0.5) +
		geom_segment(aes(x=medX, xend=medX, y=YVAR.Q1, yend=YVAR.Q3, color=CLASS), size=3, alpha=0.5) +
		geom_point(size=5,color="black", aes(x=medX, y=medY)) + geom_point(size=4, aes(x=medX, y=medY, color=CLASS)) +
		geom_label_repel(
			aes(x=medX, y=medY, fill=CLASS, label=CLASS), 
			fontface='bold',color='white',box.padding = unit(0.3, "lines"),point.padding=unit(0.6, "lines"),segment.color='grey50'
		) + theme_classic(base_size = 16) + theme(legend.position="none") + ylab(ylab) + xlab(xlab)

	if (return) return(p) else print(p)
}

.env$waterfall <- function(value, group, direction="increasing", return=FALSE) {
	require(ggplot2)
	
	df <- data.frame(Value=value, Type=group, Rank=rank(value, ties.method="random"))
	if (direction %like% "dec") df$Rank = max(df$Rank) + 1 - df$Rank
	p <- ggplot(df, aes(x=Rank, y=Value, fill=Type)) + geom_bar(stat="identity")
	if (return == TRUE) return(p) else print(p)
}

.env$mx_zscore2 <- function(x) {
 require(dplyr)
    denom <- ifelse(sd(x,na.rm=TRUE) >= 0.1,sd(x,na.rm=TRUE),0.1)
    result <- (x-mean(x,na.rm=TRUE)) / denom
    
}


#### Custom colour scheme for RAPT

theme_rapt <- function() {
    base_size <- 16
    theme_grey(
        base_size = base_size, base_family = "Source Sans Pro",
        base_line_size = base_size / 22, base_rect_size = base_size / 22
      ) %+replace%
        theme(
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid = element_line(colour = "grey95"),
          panel.grid.minor = element_line(size = rel(0.5)),
          strip.background = element_rect(fill = "grey90", colour = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=base_size, hjust = 0, family="Source Serif Pro", margin=unit(c(0,0,1,0),"lines")),
          axis.title.x = element_text(face="bold", margin=unit(c(0.75,0,0,0),"lines")),
          axis.title.y = element_text(face="bold", margin=unit(c(0,0.75,0,0),"lines"), angle=90),
          panel.background = element_rect(fill = "white", colour = NA),
          
          complete = TRUE
        )
}



rapt.extended <- c(blue="#0A95A7", orange="#E5541B", gray="#5B666F", yellow="#FE9F33", green="#4DA167", lilac="#947EB0", `dark blue`="#083D77")
rapt.darker.extended <- c(blue2="#06707E",orange2="#AC3F12",gray2="#444D54",yellow2="#D07920", green2="#356E46", lilac2="#69597D", `dark blue2`="#052344")

rapt_colors <- function(...) {
    cols <- c(...)
    full <- c(rapt.extended, rapt.darker.extended)
    if (is.null(cols)) return(full)
    return(full[cols])
}
rapt_pal <- function(palette = "main", reverse = FALSE, ...) {
    if (stringr::str_detect(palette, "dark")) {
        pal <- rapt.darker.extended
    } else if (stringr::str_detect(palette, "both")) {
        pal <-  c(rapt.extended, rapt.darker.extended)
    } else {
        pal <- rapt.extended
    }
    if (reverse) pal <- rev(pal)
    
    function(n) {
        if (n > length(pal)) stop("Palette ",palette," has a total of ",length(pal)," colors, but you need ",n,"!")
        return(unname(pal)[seq_len(n)])
    }
}

scale_color_rapt <- function(palette="main", reverse = FALSE, ...) {
    pal <- rapt_pal(palette = palette, reverse = reverse)
    discrete_scale("colour", paste0("rapt_", palette), palette = pal, ...)
}
scale_colour_rapt <- scale_color_rapt


scale_fill_rapt <- function(palette = "main", reverse = FALSE, ...) {
    pal <- rapt_pal(palette = palette, reverse = reverse)
    discrete_scale("fill", paste0("rapt_", palette), palette = pal, ...)
}

sum_mx <- function(x) {sum(x,na.rm=TRUE)}

mean_mx <- function(x) {mean(x,na.rm=TRUE)}

attach(.env)
