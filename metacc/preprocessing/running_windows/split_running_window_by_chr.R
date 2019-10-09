library(data.table)
library(purrr)

feature_metadata <- fread("/hps/nobackup/stegle/users/ricard/gastrulation/features/filt/window2000_step1000.bed", sep="\t", header=F, select=c(1,5), verbose=F, stringsAsFactors=T) %>%
  setnames(c("chr","id"))
print("foo")

data <- fread("/hps/nobackup/stegle/users/ricard/gastrulation/met/parsed/window2000_step1000.tsv", select=c(1,2,4,5), sep="\t", header=F, verbose=F, stringsAsFactors=T) %>%
	setnames(c("id_met", "id", "met_reads", "total_reads")) 
print("bar")

data <- merge(data, feature_metadata, by=c("id"))
print("baz")

data %>% split(.$chr) %>% walk2(., names(.), 
  function(x,y) fwrite(x[,c(1,2,3,4)], sprintf("/hps/nobackup/stegle/users/ricard/gastrulation/met/parsed/window_split/%s.tsv",y), quote=F, col.names=T, sep="\t" )
)

system("pigz -p 3 /hps/nobackup/stegle/users/ricard/gastrulation/met/parsed/window_split/*.tsv")


