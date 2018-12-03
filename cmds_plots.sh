

mainDirs=( "CMP_TADs_GENES_GO_pvalSelect_top" "CMP_TADs_GENES_GO_top" )


all_plots=( "GOlog10_pval" "GOfoldEnrichment" "GOgeneRatio" "g_density" "g_diameter" "g_meanEccentricity" "g_meanBetweenness" "g_meanDist" "g_meanDistDir" "intersectRatio" "signifGOidRatio" "signifGOid" "signifGOminRatio" "signifGOmin" "signifGOtermRatio" "signifGOterm" )

all_ranks=( "FCC" "avg" )

#for rank in ${all_ranks[@]}; do
#    for dir in ${mainDirs[@]}; do   
#        for toplot in ${all_plots[@]}; do

#            mkdir -p $dir/all_plots/$rank/$toplot
#            cp $dir/BP/*/$rank/*${toplot}_boxplot_nojitter.svg $dir/all_plots/$rank/$toplot
#            cp $dir/BP/*/$rank/*${toplot}_violinplot.svg $dir/all_plots/$rank/$toplot

#        done
#    done
#done



toplotvar="aucCoexprDist_aucFCC"

for rank in ${all_ranks[@]}; do
    for dir in ${mainDirs[@]}; do   
        mkdir -p $dir/all_plots/$rank/$toplotvar
        cp $dir/BP/*/$rank/*${toplotvar}.svg $dir/all_plots/$rank/$toplotvar
    done
done


#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/GO_log10
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*log10_pval_boxplot_nojitter.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots/GO_log10

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/GO_foldEnrichment
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*foldEnrichment_boxplot_nojitter.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots/GO_foldEnrichment

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/GO_geneRatio
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*geneRatio_boxplot_nojitter.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots/GO_geneRatio

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*g_density_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*g_diameter_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*g_meanEccentricity_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*g_meanBetweenness_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*g_meanDist_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*g_meanDistDir_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_intersectRatio_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_signifGOidRatio_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_signifGOid_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_signifGOminRatio_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_signifGOmin_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_signifGOtermRatio_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots

#mkdir CMP_TADs_GENES_GO_pvalSelect_top/all_plots/
#cp CMP_TADs_GENES_GO_pvalSelect_top/BP/*/*/*_signifGOterm_violinplot.svg CMP_TADs_GENES_GO_pvalSelect_top/all_plots




#mkdir CMP_TADs_GENES_GO_top/all_plots


#cp CMP_TADs_GENES_GO_top/BP/*/*/*log10_pval_boxplot_nojitter.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*foldEnrichment_boxplot_nojitter.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*geneRatio_boxplot_nojitter.svg CMP_TADs_GENES_GO_top/all_plots


#cp CMP_TADs_GENES_GO_top/BP/*/*/*g_density_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*g_diameter_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*g_meanEccentricity_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*g_meanBetweenness_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*g_meanDist_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*g_meanDistDir_violinplot.svg CMP_TADs_GENES_GO_top/all_plots


#cp CMP_TADs_GENES_GO_top/BP/*/*/*_signifGOidRatio_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*_signifGOid_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*_signifGOminRatio_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*_signifGOmin_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*_signifGOtermRatio_violinplot.svg CMP_TADs_GENES_GO_top/all_plots
#cp CMP_TADs_GENES_GO_top/BP/*/*/*_signifGOterm_violinplot.svg CMP_TADs_GENES_GO_top/all_plots


