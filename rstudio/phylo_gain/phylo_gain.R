library(caper)

calc_phylo_gain <- function(tree, tips) {
  tre <- read.tree(tree)
  clmat <- clade.matrix(tre)
  tips <- readLines(tips)
  tipsPD <- pd.calc(clmat, tip.subset = tips)
  treePD <- pd.calc(clmat)
  cat("TipsPD: ",tipsPD,'\n')
  cat("treePD: ",treePD,'\n')
  return (tipsPD/treePD)
}

print(calc_phylo_gain('mob_prots_cdhit.faa.sto.faa.tree','all_soil_plasmids_rescued_mob.txt.fixed'))

#print(calc_phylo_gain('pruned_refsoil_tree.nwk','plasmid_leaves.txt'))

#print(calc_phylo_gain('pruned_plsdb_tree.nwk','plasmid_leaves.txt'))

print(calc_phylo_gain('pruned_isolates_tree.nwk','plasmid_leaves.txt'))
