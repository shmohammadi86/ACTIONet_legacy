
map.clusters <- function(Labels, clusters) {
    N = length(Labels)
    cluster.ids = sort(unique(clusters))
    if (is.factor(Labels)) 
        Label.ids = levels(Labels) else Label.ids = sort(unique(Labels))
    
    W = sapply(Label.ids, function(label) {
        idx1 = which(Labels == label)
        n1 = length(idx1)
        if (n1 == 0) 
            return(array(0, length(cluster.ids)))
        
        log.pvals = sapply(cluster.ids, function(cluster.id) {
            idx2 = which(clusters == cluster.id)
            n2 = length(idx2)
            if (n2 == 0) 
                return(0)
            
            success = intersect(idx1, idx2)
            
            pval = phyper(length(success) - 1, n1, N - n1, n2, lower.tail = FALSE)
            return(-log10(pval))
        })
        return(log.pvals)
    })
    
    W[is.na(W)] = 0
    W[W > 300] = 300
    
    W.matched = MWM(W)
    W.matched = as(W.matched, "dgTMatrix")
    
    updated.Labels = rep(NA, length(Labels))
    matched.clusters = W.matched@i + 1
    matched.celltype = Label.ids[W.matched@j + 1]
    for (k in 1:length(matched.clusters)) {
        updated.Labels[clusters == matched.clusters[k]] = matched.celltype[[k]]
    }
    
    return(updated.Labels)
}


smooth.archetype.footprint <- function(ACTIONet.out, alpha_val = 0.9, thread_no = 8) {
    H.norm = as(t(ACTIONet.out$reconstruct.out$H_stacked), "dgTMatrix")
    cs = Matrix::colSums(H.norm)
    cs[cs == 0] = 1
    H.norm = as.matrix(sparseMatrix(i = H.norm@i + 1, j = H.norm@j + 1, x = H.norm@x/cs[H.norm@j + 1], dims = dim(H.norm)))
    
    pr.scores = batchPR(G = ACTIONet.out$build.out$ACTIONet, H.norm, alpha = alpha_val, thread_no = thread_no)
    
    arch.projections = apply(pr.scores, 2, function(x) {
        cond = sweepcut(ACTIONet.out$build.out$ACTIONet, x)
        idx = which.min(cond)
        
        perm = order(x, decreasing = TRUE)
        x[perm[(idx + 1):length(x)]] = 0
        
        return(x)
    })
    
    return(arch.projections)
}

impute.geneset.activity <- function(ACTIONet.out, sce, genes, alpha_val = 0.9, thread_no = 8) {
    imputed.gene.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, genes, alpha_val, thread_no, prune = FALSE, rescale = FALSE)
    
    
    gs.score = Matrix::rowMeans(imputed.gene.expression)
    
    return(gs.score)
}





orthoProject <- function(A, S) {
    A = scale(A)
    S = scale(S)
    A_r = A - S %*% MASS::ginv(t(S) %*% S) %*% (t(S) %*% A)
    A_r = scale(A_r)
    return(A_r)
}

is.sparseMatrix <- function(aa) {
    return(length(which(is(aa) == "sparseMatrix")) != 0)
}

# HGT tail bound
Kappa <- function(p, q) {
    kl = array(1, length(p))
    
    suppressWarnings({
        a = p * log(p/q)
    })
    a[p == 0] = 0
    
    suppressWarnings({
        b = (1 - p) * log((1 - p)/(1 - q))
    })
    b[p == 1] = 0
    
    k = a + b
    return(k)
}

HGT_tail <- function(population.size, success.count, sample.size, observed.success) {
    if (sum(success.count) == 0) 
        return(1)
    
    success.rate = success.count/population.size
    expected.success = sample.size * success.rate
    delta = (observed.success/expected.success) - 1
    
    log.tail_bound = sample.size * Kappa((1 + delta) * success.rate, success.rate)
    log.tail_bound[delta < 0] = 0
    log.tail_bound[is.na(log.tail_bound)] = 0
    
    return(log.tail_bound)
}



arch.init.clusters <- function(ACTIONet.out, update = TRUE, thread_no = 8) {
    # U = ACTIONet.out$reconstruct.out$C_stacked[, ACTIONet.out$core.out$core.archs]
    U = t(ACTIONet.out$reconstruct.out$H_stacked[ACTIONet.out$core.out$core.archs, ])
    U.smoothed = (batchPR(ACTIONet.out$build.out$ACTIONet, U, thread_no = thread_no))
    # U.smoothed = apply(U.smoothed, 2, function(x) x / max(x))
    
    initial.clusters = apply(U.smoothed, 1, which.max)
    
    if (update) {
        initial.clusters.updated = update.cell.labels(ACTIONet.out, initial.clusters)
    } else {
        initial.clusters.updated = initial.clusters
    }
    initial.clusters.updated = match(initial.clusters.updated, sort(unique(initial.clusters.updated)))
    
    return(initial.clusters.updated)
}

construct.archetype.signature.profile <- function(sce, ACTIONet.out) {
    require(SingleCellExperiment)
    
    # eigengene x archetypes
    reduced.archetype.profile = (t(sce@reducedDims[["S_r"]]) %*% ACTIONet.out$reconstruct.out$C_stacked)
    
    X = rowData(sce)
    cnames = colnames(X)
    idx = grep("^PC", cnames)
    V = as.matrix(X[, idx])
    
    perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
    V = V[, perm]
    
    # gene x archetypes
    archetype.signature.profile = V %*% reduced.archetype.profile
    rownames(archetype.signature.profile) = rownames(sce)
    
    return(archetype.signature.profile)
}

add.archetype.labels <- function(ACTIONet.out, k_min = 2, k_max = 20) {
    temp <- lapply(k_min:k_max, function(i) 1:i)
    arch.labels <- paste("A", unlist(lapply(1:length(temp), function(i) paste(i, temp[[i]], sep = "_"))), sep = "")
    
    arch.labels <- arch.labels[ACTIONet.out$reconstruct.out$selected_archs]
    colnames(ACTIONet.out$reconstruct.out$C_stacked) <- arch.labels
    rownames(ACTIONet.out$reconstruct.out$H_stacked) <- arch.labels
    colnames(ACTIONet.out$reconstruct.out$archetype_profile) <- arch.labels
    colnames(ACTIONet.out$reconstruct.out$backbone) <- arch.labels
    rownames(ACTIONet.out$reconstruct.out$backbone) <- arch.labels
    colnames(ACTIONet.out$signature.profile) <- arch.labels
    rownames(ACTIONet.out$reconstruct.out$landmark_cells) <- arch.labels
    colnames(ACTIONet.out$archetype.differential.signature) <- arch.labels
    
    return(ACTIONet.out)
}


# define utility function to adjust fill-opacity using css
fillOpacity <- function(., alpha = 0.5) {
    css <- sprintf("<style> .js-fill { fill-opacity: %s !important; } </style>", alpha)
    prependContent(., HTML(css))
}

mycircle <- function(coords, v = NULL, params) {
    library(igraph)
    
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
    }
    
    mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color, vertex.size, vertex.frame.width, FUN = function(x, y, bg, fg, 
        size, lwd) {
        symbols(x = x, y = y, bg = bg, fg = fg, lwd = lwd, circles = size, add = TRUE, inches = FALSE)
    })
}

#' @name hashid_defaults
#'
#' @title Default Values for hashid settings
#' @description Default alphabet, separators, and 
#'    ratio of character separators and guards for hashid
#'   
#' @source 
#' http://www.hashids.org
#' 

#' @rdname hashid_defaults
#' @export
DEFAULT_ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

#' @rdname hashid_defaults
#' @export
DEFAULT_SEPS = "cfhistuCFHISTU"

#' @rdname hashid_defaults
#' @export
RATIO_SEPARATORS = 3.5

#' @rdname hashid_defaults
#' @export
RATIO_GUARDS = 12

#' Decodes a hashid into the original integer or integer vector
#'
#' @param hash_str hashid string to decode into integer or integer vector
#' @param settings Settings list generated by hashid_settings
#' @return integer or integer vector
#' @export
#'
DECODE = function(hash_str, settings) {
    if (hash_str == "") 
        stop("decode: invalid hashid")
    
    salt = settings$salt
    alphabet = settings$alphabet
    separator = settings$separator
    guards = settings$guards
    
    parts = SPLIT(hash_str, guards)
    hashid = ifelse(2 <= length(parts) & length(parts) <= 3, parts[2], parts[1])
    
    if (hashid == "") 
        stop("decode: invalid hashid, cannot decode")
    lottery = substr(hashid, 1, 1)
    hashid = substr(hashid, 2, nchar(hashid))
    
    hash_parts = SPLIT(hashid, separator)
    unhashed_parts = c()
    for (p in hash_parts) {
        alphabet_salt = substr(paste0(lottery, salt, alphabet), 1, nchar(alphabet))
        alphabet = shuffle(alphabet, alphabet_salt)
        unhashed_parts = c(unhashed_parts, unhash(p, alphabet))
    }
    
    rehash = tryCatch({
        ENCODE(unhashed_parts, settings)
    }, error = function(e) {
        stop("decode: invalid hashid, cannot decode")
    })
    if (!all(hash_str == rehash)) {
        stop("decode: invalid hashid, cannot decode")
    }
    
    return(unhashed_parts)
}

#' Encodes an integer or integer vector into a hashid string.
#' All numbers must be non-negative integers.
#'
#' @param int Integer or integer vector to encode
#' @param settings Settings list generated by hashid_settings
#' @return hashid string
#' @export
#'
ENCODE = function(int, settings) {
    if (!all(c("alphabet", "salt", "guards", "separator", "min_length") %in% names(settings))) {
        stop("encode: missing some parameters in settings list")
    }
    if (any(int < 0)) {
        stop("encode: numbers must be non-negative")
    }
    if (any(int%%1 != 0)) {
        stop("encode: numbers must be integers")
    }
    if (length(int) < 1) {
        stop("encode: Invalid length!")
    }
    
    alphabet = settings$alphabet
    salt = settings$salt
    guards = settings$guards
    separator = settings$separator
    min_length = settings$min_length
    alphabet_len = nchar(settings$alphabet)
    sep_len = nchar(settings$separator)
    
    vec_hash = sum(sapply(1:length(int), function(i) {
        int[i]%%(100 + i - 1)
    }))
    # lottery character
    lottery = substr(alphabet, (vec_hash%%alphabet_len) + 1, (vec_hash%%alphabet_len) + 1)
    encoded = lottery
    
    for (i in 1:length(int)) {
        alphabet_salt = substr(paste0(lottery, salt, alphabet), 1, alphabet_len)
        alphabet = shuffle(alphabet, alphabet_salt)
        last = hash(int[i], alphabet)
        encoded = paste0(encoded, last)
        int[i] = int[i]%%(ascii_val(substr(last, 1, 1)) + (i - 1))
        encoded = paste0(encoded, substr(separator, (int[i]%%sep_len + 1), (int[i]%%sep_len + 1)))
    }
    
    encoded = substr(encoded, 1, nchar(encoded) - 1)
    if (nchar(encoded) <= min_length) {
        encoded = enforce_min_length(encoded, min_length, alphabet, guards, vec_hash)
    }
    
    return(encoded)
}

#' A function to create a hashid settings list.  
#'
#' @param salt An additional string to make hashids more unique.
#' @param min_length Minimum length for hashid.
#' @param alphabet String of characters for hashid.
#' @param sep String of characters to use as separators.
#' @return A list of parameters used in encoding and decoding.
#' @export
#'
hashid_settings = function(salt, min_length = 0, alphabet = DEFAULT_ALPHABET, sep = DEFAULT_SEPS) {
    
    alphabet_vec = unique(strsplit(alphabet, split = "")[[1]])
    sep_vec = unique(strsplit(sep, split = "")[[1]])
    
    separator_ = paste(intersect(sep_vec, alphabet_vec), collapse = "")
    alphabet_ = paste(setdiff(alphabet_vec, sep_vec), collapse = "")
    
    if (nchar(separator_) + nchar(alphabet_) < 16) {
        # if(nchar(alphabet_) < 16) {
        stop("hashid_settings: Alphabet must be at least 16 unique characters.")
    }
    
    separator_ = shuffle(separator_, salt)
    min_separators = ceiling(nchar(alphabet_)/RATIO_SEPARATORS)
    
    ## if needed get more separators from alphabet ##
    if (nchar(separator_) < min_separators) {
        if (min_separators == 1) 
            min_separators = 2
        split_at = min_separators - nchar(separator_)
        separator_ = paste0(separator_, substr(alphabet_, 1, split_at))
        alphabet_ = substr(alphabet_, split_at + 1, nchar(alphabet_))
    }
    
    alphabet_ = shuffle(alphabet_, salt)
    num_guards = ceiling(nchar(alphabet_)/RATIO_GUARDS)
    
    if (nchar(alphabet_) < 3) {
        guards_ = substring(separator_, 1, num_guards)
        separator_ = substr(separator_, num_guards + 1, nchar(separator_))
    } else {
        guards_ = substring(alphabet_, 1, num_guards)
        alphabet_ = substr(alphabet_, num_guards + 1, nchar(alphabet_))
    }
    
    return(list(alphabet = alphabet_, salt = salt, guards = guards_, separator = separator_, min_length = min_length))
}

#' Calculate the ascii value number of a character
#' 
#' @param char character
#' @return ascii value integer
#' 
ascii_val = function(char) {
    if (!is.character(char)) 
        stop("ascii_val: must be character")
    strtoi(charToRaw(char), 16)
}

#' Converts a base 16 string to a base 10 number.
#' Because I couldn't get base R functions to work for big hex numbers.
#' 
#' @param str_16 base 16 number as a string.
#' @return base 10 integer.
#' 
base16_to_dec = function(str_16) {
    str_vec = strsplit(tolower(str_16), split = "")[[1]]
    str_vec = sapply(str_vec, function(x) {
        if (x %in% as.character(0:9)) {
            as.numeric(x)
        } else if (x %in% c("a", "b", "c", "d", "e", "f")) {
            ascii_val(x) - 87
        } else {
            stop("base16_to_dec: Invalid hex character")
        }
    })
    
    vec_pwrs = 16^(rev(1:length(str_vec)) - 1)
    
    sum(vec_pwrs * str_vec)
}

#' Converts a base 10 number to base 16 number.  
#' Because I couldn't get R's as.hexmode() to work for big integers.
#'
#' @param dec base 10 integer
#' @return base 16 number as a string
#'
dec_to_base16 = function(dec) {
    num_vec = c()
    while (dec > 0) {
        rem = dec%%16
        num_vec = c(rem, num_vec)
        dec = floor(dec/16)
    }
    
    hex_vec = sapply(num_vec, function(x) {
        if (x < 10) {
            return(x)
        } else {
            base::letters[x - 9]
        }
    })
    
    paste(hex_vec, collapse = "")
}

#' Enforces hashid minimum length by padding the hashid with additional characters.
#'
#' @param encoded encoded hashid
#' @param min_length minimum length required for hashid
#' @param alphabet set of letters used to generate hashid
#' @param guards set of guards used to generate hashid
#' @param values_hash value hashed used to select guard characters
#' @return hashid with padded characters to insure minimum length
#'
enforce_min_length = function(encoded, min_length, alphabet, guards, values_hash) {
    
    guards_len = nchar(guards)
    guards_idx = (values_hash + ascii_val(substr(encoded, 1, 1)))%%guards_len + 1
    encoded = paste0(substr(guards, guards_idx, guards_idx), encoded)
    
    if (nchar(encoded) < min_length) {
        guards_idx = (values_hash + ascii_val(substr(encoded, 3, 3)))%%guards_len + 1
        encoded = paste0(encoded, substr(guards, guards_idx, guards_idx))
    }
    
    split_at = nchar(alphabet)/2 + 1
    while (nchar(encoded) < min_length) {
        alphabet = shuffle(alphabet, alphabet)
        encoded = paste0(substr(alphabet, split_at, nchar(alphabet)), encoded, substr(alphabet, 1, split_at - 1))
        excess = nchar(encoded) - min_length
        
        if (excess > 0) {
            from_index = floor(excess/2) + 1
            encoded = substr(encoded, from_index, from_index + min_length - 1)
        }
    }
    
    return(encoded)
}

#' Maps an integer to a string.
#' Generated string will be inversely proportional to alphabet length.
#' 
#' @param number Integer to hash
#' @param alphabet Possible letters for string.
#' @return hashed string
#'
hash = function(number, alphabet) {
    alphabet_len = nchar(alphabet)
    alphabet_vec = strsplit(alphabet, split = "")[[1]]
    
    hashed = c()
    while (number > 0) {
        hash_idx = (number%%alphabet_len) + 1
        hashed = c(alphabet_vec[hash_idx], hashed)
        number = floor(number/alphabet_len)
    }
    
    return(paste(hashed, collapse = ""))
}

#' Permutes the characters in a string based on an inputted salt string.
#'
#' @param string String to be permuted
#' @param salt cryptograph salt string that is used to permute strings
#' @return shuffled string
#'
shuffle = function(string, salt) {
    salt_len = nchar(salt)
    str_len = nchar(string)
    if (salt_len < 1 | str_len < 2) 
        return(string)
    
    salt_sum = 0
    salt_index = 1
    string_vec = strsplit(string, split = "")[[1]]
    salt_vec = strsplit(salt, split = "")[[1]]
    for (i_str in rev(2:str_len)) {
        ## Pseudo Randomize based on salt ##
        salt_index = (salt_index - 1)%%(salt_len) + 1
        salt_int_val = ascii_val(salt_vec[salt_index])
        salt_sum = salt_sum + salt_int_val
        swap_pos = (salt_sum + salt_index + salt_int_val - 1)%%(i_str - 1) + 1
        
        ## Swap positions ##
        temp = string_vec[swap_pos]
        string_vec[swap_pos] = string_vec[i_str]
        string_vec[i_str] = temp
        salt_index = salt_index + 1
    }
    
    return(paste(string_vec, collapse = ""))
}

#' Splits a string based on a set of splitting characters
#'
#' @param string String to split
#' @param splitters set of splitting characters as a string
#' @return split vector of characters
#' 
SPLIT = function(string, splitters) {
    string_vec = strsplit(string, split = "")[[1]]
    split_vec = strsplit(splitters, split = "")[[1]]
    
    word = ""
    words = c()
    for (i in 1:length(string_vec)) {
        if (string_vec[i] %in% split_vec) {
            words = c(words, word)
            word = ""
        } else {
            word = paste0(word, string_vec[i])
        }
    }
    words = c(words, word)
    
    return(words)
}

#' Unhashes a string to an integer based on alphabet.
#' 
#' @param hashed String to unhash
#' @param alphabet Set of letters used for hashing
#' @return Unhashed integer
#'
unhash = function(hashed, alphabet) {
    hashed_len = nchar(hashed)
    alphabet_len = nchar(alphabet)
    alphabet_vec = strsplit(alphabet, split = "")[[1]]
    hashed_vec = strsplit(hashed, split = "")[[1]]
    
    number = 0
    for (i in 1:hashed_len) {
        position = which(alphabet_vec == hashed_vec[i]) - 1
        number = number + (position * alphabet_len^(hashed_len - i))
    }
    
    return(number)
}

textHalo <- function(x, y = NULL, labels, col = "white", bg = "black", r = 0.1, ...) {
    
    theta = seq(0, 2 * pi, length.out = 50)
    xy <- xy.coords(x, y)
    xo <- r * strwidth("A")
    yo <- r * strheight("A")
    
    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text(xy$x + cos(i) * xo, xy$y + sin(i) * yo, labels, col = bg, ...)
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col = col, ...)
}


combine.logPvals <- function(logPvals, top.len = NULL, base = 10) {
    if (is.null(top.len)) {
        top.len = nrow(logPvals)
    }
    kappa = 1/log(exp(1), base = base)
    logPvals = kappa * logPvals
    
    combbined.log.pvals = -apply(logPvals, 2, function(lx) {
        perm = order(lx, decreasing = T)
        
        return(log(top.len) - logSumExp(lx[perm[1:top.len]]))
    })
    
    return(combbined.log.pvals)
}

