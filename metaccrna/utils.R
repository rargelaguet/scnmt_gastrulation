# return id's from data.table sorted based on lci (decreasing)
lci_sorted_feats <- function(dt=met_dt_list$Promoters) {
  dt %>% .[,c("id", "anno", "lci")] %>% 
    .[!duplicated(.)] %>% 
    .[order(-lci)] %>% 
    .$id
}

weight_se <- function(r, total ) return( 1/(r*(1-r)/(total)) )

calc_site_stats <- function(filePath=NULL,
                            keep_samples = NULL, ## pass QC samples
                            min_N = 3, ## min number of calls at site
                            min_cov = 500, ## min number of cells detecting
                            alpha = 0.1) ## for lower bound of CI
{
  dt <- fread(filePath, colClasses=c("factor","character","factor","integer","integer","numeric")) %>% 
    setnames(c("sample","id","anno","Nmet","N","rate")) 
  
  if (!is.null(keep_samples))
    dt <-  dt[sample %in% keep_samples]
  
  dt %>% 
    # Filter by coverage
    .[N >= min_N] %>% # within a region
    .[,cell_cov:=.N, by=c("anno", "id")] %>% .[cell_cov >= min_cov] %>% .[,cell_cov:=NULL] %>% # across cells
    .[,rate:=(Nmet+1)/(N+2)] %>% # MAP estimate for ratex
    .[,wij:=weight_se(rate, N)] %>%  # weights based on SE of MAP estimates
    .[,rbar:=sum(wij*rate)/sum(wij), by=c('anno', 'id')] %>%  # mean across cells
    .[,n_i:= (sum(wij)^2 - sum(wij^2))/sum(wij), by=c('anno', 'id')] %>%  # sum of sample weights at region
    .[,vhat:= sum(wij*(rate-rbar)^2)/n_i,  by=c('anno', 'id')] %>%  # site variance across cells
    .[,lci:= n_i*vhat/(qchisq(p=1-alpha/2, df = n_i)),  by=c('id', 'anno')] # lower bound of CI
}


perm <- function(cell_names=c("foo", "bar", "egg", "spam", "xyz"), ids = paste0("id",1:2)) {
  lapply(seq_len(length(cell_names)-1), function(i) {
    ## ----------- Permutation of cells and ids
    ## helper function to create unique permutations of cells and ids for MDS distances:
    
    expand.grid(cell1=cell_names[i], cell2=cell_names[-1:-i], id = ids, stringsAsFactors = FALSE)
  }) %>% data.table::rbindlist() 
}


wt_euc_dist <- function(met_dt = met_dt, 
                        sample_name = "id_met",
                        n_hvr = NULL ## number of highly variable regions to use
                        ## , subset=NULL ## seems to throw error
) {
  ## ----------- Weighted Euclidean Distance
  ## 
  ## methods: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770512/
  if (!is.null(n_hvr)) {
    all_ids <- lci_sorted_feats(met_dt)
    keep_ids <- all_ids[seq_len(min(n_hvr, length(all_ids)))]
    met_dt <- met_dt[id %in% keep_ids]
  }
  ## pairwise weights for all cells for each id
  all_cells <- unique(met_dt[,sample_name, with=FALSE]) %>% unlist() %>% sort()
  iter <- perm(cell_names = all_cells, ids = unique(met_dt$id))
  ## to use for joins
  joiner <- met_dt[,c(sample_name, "id", "rate", "wij", "N"), with=FALSE] %>% set_names(c("cell1", "id", "wij", "rate", "N"))
  
  ## add cell1 weight and N
  iter <- merge(iter, joiner, by = c("id", "cell1"), all.x = TRUE)
  
  ## add cell2 weight and N
  iter <- merge(iter, joiner %>% setnames("cell1", "cell2"), by = c("id", "cell2"), all.x = TRUE)
  
  ## add site total bi-cell weight as sqrt of products, if both present
  iter[,wi_jk:=sqrt(wij.x*wij.y)]
  
  ## sum of weights
  iter[,sumwt:=sum(wi_jk, na.rm = TRUE), by=id]
  iter[,sumN:=sum(N.x+N.y, na.rm = TRUE), by=id]
  ## make some of the weights to add u to total observations
  ## this probably can be percentage of CpGs covered
  ## although looking at single id's there's no difference, just wondering
  ## if it implies more observations for an id makes that id more important? ##TODO
  iter[,wi_jk:=wi_jk*sumN/sumwt]
  
  ## weighted euc dist
  iter[,l2_wt:=wi_jk*(rate.x-rate.y)^2]
  ## pairwise distances
  dists <- iter[, .(wt_euc_d=sqrt(sum(l2_wt, na.rm = TRUE))), by=c('cell1', 'cell2')]
  
  ## add same-cell rows as well for dcast to create dist mat easily
  dists <- rbind(dists, data.table(cell1 = all_cells, cell2=all_cells, wt_euc_d=0))
  ## get the distance matrix
  dist_mat  <- dcast(dists, cell1~cell2, value.var = "wt_euc_d") %>%
    data.frame(row.names = TRUE, check.names = FALSE)
  ## order dim names
  dist_mat <- dist_mat[all_cells, all_cells]
  stopifnot(all(is.na(dist_mat[lower.tri(dist_mat)])))
  ## copy upper tri to lower tri
  dist_mat[lower.tri(dist_mat)] <- 0
  dist_mat <- dist_mat+t(dist_mat)
  
  
  return(dist_mat)
}


get_dist_mats <- function(dt_list, cores=parallel::detectCores(), 
                          sample_name = "id_met",
                          n_hvr = 100
                          ## , subset=NULL ## seems to throw error
) {
  ## ----------- Calculate Distance Matrices 
  
  require(parallel)
  ## check that they're all data.tables
  stopifnot(all(sapply(dt_list, function(x) is(x, "data.table"))))
  
  mclapply(dt_list, function(x) {
    rt <- system.time({
      dm <- wt_euc_dist(met_dt = x, 
                        sample_name = sample_name,
                        n_hvr = n_hvr)
    })["elapsed"]
    
    return(list(runtime=rt, distmat=dm))
  }, mc.cores = cores)
}


## ----------- dcast functions ----------- 
## sort dt that has lci, anno, id columns based on lci and retun sorted id's in decreasing order

## dcast 'id' against sample_col for value_var, make 'id' rowname, sort by lci, create matrix
dcast_sort <- function(w, 
                       value_var="rate", 
                       sample_col = "id_met") {
  dcast.data.table(w,  as.formula(sprintf("id~%s", sample_col)), 
                   value.var = value_var) %>% 
    data.frame(row.names = TRUE) %>%
    .[lci_sorted_feats(w),] %>% as.matrix()
}

## no sorting
dcast_by <- function(w=met_dt, 
                     value_var="rate", 
                     sample_col = "id_met") {
  dcast(w,  as.formula(sprintf("id~%s", sample_col)), 
        value.var = value_var) %>% 
    data.frame(row.names = TRUE, check.names = FALSE)
}


dcast_dt_list <-
  function(dt_list,
           row = "id",
           col = "id_met",
           value.var = "rate") {
    ## ----------- dcast list of data.tables to list of matrices 
    lapply(dt_list,
           function(w) {
             dcast.data.table(w, as.formula(sprintf("%s~%s", row, col)), value.var = value.var) %>%
               data.frame(row.names = TRUE, check.names = FALSE) %>%
               .[lci_sorted_feats(w),] %>%
               as.matrix()
           })
  }

## ----------- rbindListWithNames ----------- 
rbindListWithNames <- function(lst, new_col = "dataset") {
  lst_with_newcol <- mapply(x=names(lst), y=lst, FUN = function(x, y){
    y[,new_col] <- x
    y
  }, SIMPLIFY = FALSE)
  Reduce(rbind, lst_with_newcol)
}

## ----------- installer ----------- 
installer <- function(pkgs) {
  sapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat('installing: ', pkg, '\n')
      BiocManager::install(pkg)
    }
    NULL
  })
  NULL
}