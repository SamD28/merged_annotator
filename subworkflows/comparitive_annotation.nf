#!/usr/bin/env nextflow

process create_metadata_summary {
    publishDir "${params.outdir}/annotation_results", mode: 'copy'

    input:
    path metadata_file
    output:
    path "metadata_column_COMPARATIVE_ANNOTATION_summary.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv("${metadata_file}", sep=',', header=0)

    summary = pd.DataFrame({
        'Column_Index': range(len(df.columns)),
        'Column_Name': df.columns
    })

    summary.to_csv('metadata_column_COMPARATIVE_ANNOTATION_summary.tsv', sep='\t', index=False)
    """
}

process run_skani_annotation {
    publishDir "${params.annotation_results}/ani", mode: 'copy'
    cpus params.cpus

    input:
    path(genome_dir)

    output:
    path "skani_fullmatrix"
    path "skani_ANI_dist.tsv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    skani triangle *.fa -t ${task.cpus} --full-matrix | tail -n +2 | sed -e 's/.fa//g' > skani_fullmatrix

    skani triangle *.fa -t ${task.cpus} --full-matrix --distance > test
    #skani triangle *.fa -t ${task.cpus} --full-matrix -l rl --distance > test
    cut -f1 test | tail -n+2 | tr '\n' '\t' | sed -e 's/^/\t/g' | sed -e 's/\t\$/\\n/' > skani_ANI_dist.tsv
    tail -n+2 test >> skani_ANI_dist.tsv
    sed -i 's/\\.fa//g' skani_ANI_dist.tsv
    rm test


    """
}

process run_skani_visualization {
    publishDir "${params.visualization_results}/ani/column_${params.metacol}", mode: 'copy'

    input:
    path(skani_fullmatrix)
    path(metadata_file)

    output:
    path "heatmap_skani.pdf"
    path "skani_interactive"

    when: params.run_visualization

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/skani_visualization.R \
    -i ${skani_fullmatrix} -m ${metadata_file} -mc ${params.metacol} -out heatmap_skani.pdf -out_html skani_interactive
    """
}

process run_skani {
    publishDir "${params.outdir}/ani" , mode : 'copy'

    cpus params.cpus

    input:
    path(genome_dir)
    path metadata_file 
    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    skani triangle *.fa -t ${task.cpus}  --full-matrix | tail -n +2 | sed -e 's/.fa//g' > skani_fullmatrix
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/skani_visualization.R \
    -i skani_fullmatrix -m ${metadata_file}  -mc ${params.metacol} -out heatmap_skani.pdf -out_html skani_interactive

    """
}

process run_panaroo {
    publishDir "${params.outdir}" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"

    input:
    path(prokka_dir)
    output:
    path "*"

    script:
    """
    source activate COMPARATIVE_ANNOTATION
    panaroo -i */*.gff -o panaroo_result --clean-mode moderate --remove-invalid-genes \
    -c 0.9 -f 0.5 --merge_paralogs --threads ${task.cpus}
    """
}


process run_ppanggolin {
    publishDir "${params.outdir}" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(prokka_dir)
    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    ls */*.gff | while read gff_file; do 
    printf "%s\t%s\n" "\$(basename "\$gff_file" .gff)" "\$(realpath "\$gff_file")"
    done > gff_paths.tsv
    ppanggolin 
    ppanggolin all --anno gff_paths.tsv -o ppanggolin_result --cpu ${task.cpus}\
    --coverage ${params.pan_identity} --identity ${params.pan_coverage} \
    --kingdom  ${params.kingdom} --rarefaction
    cd ppanggolin_result 
    ppanggolin fasta -p pangenome.h5 --prot_families all -o pan_genes  --genes all 
    python ${params.scripts_baseDir}/make_gene_count_table_ppanggolin205.py ./
    # genome_id_linking.tsv, gene_count_matrix.tsv is generated 
    """
}
//#panaroo -i */*.gff -o ppanggorin_result --clean-mode moderate --remove-invalid-genes \
//    #-c 0.9 -f 0.5 --merge_paralogs --threads ${params.cpus}
//find "\$(pwd)" -type f -name "*.gff" -exec sh -c 'echo -e "\$(basename {} .gff)\t\$(realpath "{}")"' \; > gff_paths.tsv    

process run_genePA_cluster {
    publishDir "${params.outdir}/genePA_cluster" , mode : 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    path metadata 
    output:
    path("pcoa_plot_interactive.html")    

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate scoary-2
    python ${params.scripts_baseDir}/genome_PA_stat.py --cpus ${task.cpus} --metadata ${params.metadata}

    """

}


process run_kofamscan_annotation {
    publishDir "${params.annotation_results}/kofamscan", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("ko_matrix.csv")
    path("KO_definition_GeneID_countgenomes.csv")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    
    echo running kofamscan
    time exec_annotation \
    -o kofamscan_result --cpu ${task.cpus} --e-value ${params.kofamscan_eval}\
    -c ${params.db_baseDir}/kofam/27Nov2023/config-template.yml \
    -f mapper ppanggolin_result/pan_genes/all_protein_families.faa
    
    python ${params.scripts_baseDir}/kofam_to_geneID_KO_ppanggolin205.py \
    -i ppanggolin_result/matrix.csv -o output_gene_pa_ko.csv
#KO_definition_GeneID_countgenomes.py generate KO_definition_GeneID_countgenomes.csv 
    python ${params.scripts_baseDir}/KO_definition_GeneID_countgenomes.py -i  output_gene_pa_ko.csv \
    -o KO_definition_GeneID_countgenomes.csv -k ${params.db_baseDir}/kofam/27Nov2023/ko_list

    python ${params.scripts_baseDir}/KO_Genome_long2matrix_ppanggolin205.py \
    -i output_gene_pa_ko.csv -o ko_matrix.csv
    """
}

process run_kofamscan_visualization {
    publishDir "${params.visualization_results}/kofamscan", mode: 'copy'
    publishDir "${params.visualization_results}/kofamscan/column_${params.metacol}", mode: 'copy'
    
    input:
    path(ko_matrix)
    path(metadata_file)
    
    output:
    path("KEGG_module_visualization_shiny")
    path("KEGG_module_completeness.csv")
    path("heatmap_KEGG.pdf")

    when: params.run_visualization

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/KO_module_visualization_apptainer.R \
    -i ${ko_matrix} -m ${metadata_file} -n ${params.metacol} -mc ${params.module_completeness} \
    -out_table KEGG_module_completeness.csv -out_html KEGG_module_visualization_shiny
    """
}


process run_kofamscan {
    publishDir "${params.outdir}/kofamscan" , mode : 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    path metadata_file 

    output:
    path("ko_matrix.csv")
    path("KEGG_module_visualization_shiny")
    path("KEGG_module_completeness.csv")
    path("heatmap_KEGG.pdf")

    script:
    """
    echo running kofamscan
    time exec_annotation \
    -o kofamscan_result --cpu ${task.cpus}  --e-value 0.00001 \
    -c ${params.db_baseDir}/kofam/27Nov2023/config-template.yml \
    -f mapper ppanggolin_result/pan_genes/all_protein_families.faa
    
    head kofamscan_result
    

    python ${params.scripts_baseDir}/kofam_to_geneID_KO_ppanggolin205.py \
    -i ppanggolin_result/matrix.csv -o output_gene_pa_ko.csv

    python ${params.scripts_baseDir}/KO_Genome_long2matrix_ppanggolin205.py \
    -i output_gene_pa_ko.csv  -o ko_matrix.csv

    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/KO_module_visualization_apptainer.R \
    -i  ko_matrix.csv -m ${metadata_file}  -n ${params.metacol} -mc ${params.module_completeness}\
    -out_table KEGG_module_completeness.csv -out_html KEGG_module_visualization_shiny
    """
}

process run_VFDB_annotation {
    publishDir "${params.annotation_results}/VFDB", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("pangene_vfdb_result.txt"), emit: vfdb_result
    path("gene_PA_VFDB_added.csv"), emit: vfdb_gene_pa    
    path("gene_count_VFDB_added.csv"), emit: vfdb_gene_count

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    micromamba activate rgi603
    diamond blastp -d ${params.db_baseDir}/VFDB/VFDB_setB_prot_Aug2023.dmnd \
    -q ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_vfdb_result.txt -f 6 --max-target-seqs 1 --id ${params.VFDB_identity} \
    --subject-cover ${params.VFDB_coverage} -e ${params.VFDB_e_value} 


    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt -g ppanggolin_result/gene_presence_absence.Rtab \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_PA_VFDB_added.csv


    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt \
    -g ppanggolin_result/gene_count_matrix.tsv \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_count_VFDB_added.csv     
    """
}

process run_VFDB_visualization {
    publishDir "${params.visualization_results}/VFDB", mode: 'copy'

    input:
    path(gene_PA_VFDB_added)
    path(gene_count_VFDB_added)
    path(metadata_file)

    output:
    path("heatmap_VFDB_gene_PA_metadata${params.metacol}th_col.pdf")
    path("VFDB_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_VFDB_gene_count_metadata${params.metacol}th_col.pdf")
    path("VFDB_interactive_gene_count_metadata${params.metacol}th_col")

    when: params.run_visualization

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
        # Run visualization for gene_PA
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i ${gene_PA_VFDB_added} -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_VFDB.pdf heatmap_VFDB_gene_PA_metadata${params.metacol}th_col.pdf

    # Run visualization for gene_count
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i ${gene_count_VFDB_added} -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_VFDB.pdf heatmap_VFDB_gene_count_metadata${params.metacol}th_col.pdf
    """
}

process run_VFDB {
    publishDir "${params.outdir}/VFDB" , mode : 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    path metadata_file 


    output:
    path("pangene_vfdb_result.txt")
    path("gene_PA_VFDB_added.csv")
    path("heatmap_VFDB.pdf")
    path("VFDB_interactive")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    #transeq panaroo_result/pan_genome_reference.fa pan_genome_reference.faa
    #sed -i "s/_[0-9]\$//g" pan_genome_reference.faa
    #ppanggolin_result/pan_genes/all_protein_families.faa
    micromamba activate rgi603
    diamond blastp -d ${params.db_baseDir}/VFDB/VFDB_setB_prot_Aug2023.dmnd \
    -q ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_vfdb_result.txt -f 6 --max-target-seqs 1 --id 50 \
    --subject-cover 80 -e 1e-10 

    python ${params.scripts_baseDir}/add_VFDB_togenePA.py \
    -d pangene_vfdb_result.txt -g ppanggolin_result/gene_presence_absence.Rtab \
    -va ${params.db_baseDir}/VFDB/VFDB_setB_all_informatoin.tsv \
    -o gene_PA_VFDB_added.csv

    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/VFDB_visualization.R \
    -i  gene_PA_VFDB_added.csv -m ${metadata_file} -mc ${params.metacol} -out_html VFDB_interactive

    """
}

process run_rgi_CARD_annotation {
    publishDir "${params.annotation_results}/CARD", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("pangene_rgi_CARD_result.txt"), emit: rgi_result
    path("gene_PA_CARD_added.csv"), emit: card_gene_pa
    path("gene_count_CARD_added.csv"), emit: card_gene_count
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate COMPARATIVE_ANNOTATION
    micromamba activate rgi603
    rgi load -i ${params.db_baseDir}/CARD/card.json --local 
    rgi main -i ppanggolin_result/pan_genes/all_protein_families.faa -o pangene_rgi_CARD_result --clean  \
    -t protein -n ${task.cpus} --include_nudge --local 

    python ${params.scripts_baseDir}/add_rgi_togenePA.py -i pangene_rgi_CARD_result.txt \
    -o gene_PA_CARD_added.csv -r ${params.db_baseDir}/CARD/aro_index.tsv \
    -gpa ppanggolin_result/gene_presence_absence.Rtab

    python ${params.scripts_baseDir}/add_rgi_togenePA.py -i pangene_rgi_CARD_result.txt \
    -o gene_count_CARD_added.csv -r ${params.db_baseDir}/CARD/aro_index.tsv \
    -gpa ppanggolin_result/gene_count_matrix.tsv    
    """
}

process run_rgi_CARD_visualization {
    publishDir "${params.visualization_results}/CARD", mode: 'copy'

    input:
    path(gene_PA_CARD_added)
    path(gene_count_CARD_added)
    path(metadata_file)

    output:
    path("heatmap_CARD_gene_PA_metadata${params.metacol}th_col.pdf")
    path("CARD_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_CARD_gene_count_metadata${params.metacol}th_col.pdf")
    path("CARD_interactive_gene_count_metadata${params.metacol}th_col")

    when: params.run_visualization

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    
    # Run visualization for gene_PA
    Rscript ${params.scripts_baseDir}/CARD_visualization.R \
    -i ${gene_PA_CARD_added} -m ${metadata_file} -mc ${params.metacol} -out_html CARD_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_CARD.pdf heatmap_CARD_gene_PA_metadata${params.metacol}th_col.pdf

    # Run visualization for gene_count
    Rscript ${params.scripts_baseDir}/CARD_visualization.R \
    -i ${gene_count_CARD_added} -m ${metadata_file} -mc ${params.metacol} -out_html CARD_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_CARD.pdf heatmap_CARD_gene_count_metadata${params.metacol}th_col.pdf
    """
}

process run_defensefinder_annotation {
    publishDir "${params.annotation_results}/defensefinder", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("all_protein_families_defense_finder_genes.tsv"), emit: defensefinder_result
    path("gene_PA_defense.csv"), emit: defense_gene_pa
    path("gene_count_defense.csv"), emit: defense_gene_count

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate defensefinder

    defense-finder run -w ${task.cpus} --db-type unordered --models-dir ${params.db_baseDir}/defensefinder/models --antidefensefinder -o defensefinder_out --db-type unordered ppanggolin_result/pan_genes/all_protein_families.faa

    mv defensefinder_out/all_protein_families_defense_finder_genes.tsv .
    micromamba activate COMPARATIVE_ANNOTATION
    python ${params.scripts_baseDir}/gene_matrix_add_defense.py -d all_protein_families_defense_finder_genes.tsv -g ppanggolin_result/gene_presence_absence.Rtab -o gene_PA_defense.csv

    python ${params.scripts_baseDir}/gene_matrix_add_defense.py -d all_protein_families_defense_finder_genes.tsv -g ppanggolin_result/gene_count_matrix.tsv -o gene_count_defense.csv
    """
}

process run_defensefinder_visualization {
    publishDir "${params.visualization_results}/defensefinder", mode: 'copy'

    input:
    path(gene_PA_defense)
    path(gene_count_defense)
    path(metadata_file)

    output:
    path("heatmap_defensefinder_gene_PA_metadata${params.metacol}th_col.pdf")
    path("defensefinder_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_defensefinder_gene_count_metadata${params.metacol}th_col.pdf")
    path("defensefinder_interactive_gene_count_metadata${params.metacol}th_col")

    when: params.run_visualization

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    
    # Run visualization for gene_PA
    Rscript ${params.scripts_baseDir}/Defense_visualization.R \
    -i ${gene_PA_defense} -m ${metadata_file} -mc ${params.metacol} -out_html defensefinder_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_defensefinder.pdf heatmap_defensefinder_gene_PA_metadata${params.metacol}th_col.pdf

    # Run visualization for gene_count
    Rscript ${params.scripts_baseDir}/Defense_visualization.R \
    -i ${gene_count_defense} -m ${metadata_file} -mc ${params.metacol} -out_html defensefinder_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_defensefinder.pdf heatmap_defensefinder_gene_count_metadata${params.metacol}th_col.pdf
    """
}


process run_dbCAN_annotation_bak {
    publishDir "${params.annotation_results}/dbCAN", mode: 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("dbcan_hmmerfamily_count_matrix.csv")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate run_dbcan4
    export OMP_NUM_THREADS=${task.cpus}

    time run_dbcan ppanggolin_result/pan_genes/all_protein_families.faa protein --dia_cpu ${task.cpus} \
    --hmm_cpu ${task.cpus} --out_dir db_can_out --db_dir ${params.db_baseDir}/dbCAN

    cut -f1,3 db_can_out/overview.txt | grep -Pv "\\t-\$" | awk -F'\\t' '{split(\$2,a,"+"); for (i in a) {print \$1 "\\t" a[i];}}' | \
    awk -F"\\t" '{split(\$2,a,"_"); \$2=a[1]; print \$1 "\\t" \$2}' | sed -e "s/(.*//" > dbcan_family.tsv

    python ${params.scripts_baseDir}/dbcan_result_parse.py
    """
}


process run_dbCAN_annotation {
    publishDir "${params.annotation_results}/dbCAN", mode: 'copy'
    cpus params.cpus  // 

    input:
    path(ppanggolin_dir)

    output:
    path("db_can_out"), emit: db_can_out
    path("dbcan_raw_gene_count_data.csv"), emit: raw_gene_count_data
    path("dbcan_geneID_HMMER_count_gene_count.csv"), emit: geneID_HMMER_count_gene_count
    path("dbcan_HMMER_count_gene_count_matrix.csv"), emit: HMMER_count_gene_count_matrix
    path("dbcan_raw_gene_PA_data.csv"), emit: raw_gene_PA_data
    path("dbcan_geneID_HMMER_count_gene_PA.csv"), emit: geneID_HMMER_count_gene_PA
    path("dbcan_HMMER_count_gene_PA_matrix.csv"), emit: HMMER_count_gene_PA_matrix

    script:
    def parallel_jobs = (params.cpus as int).intdiv(4)  // divide maximum cpu by 4 
    """
    # Activate environment for file splitting
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC

    # Split the input file
    seqkit split2 ppanggolin_result/pan_genes/all_protein_families.faa -p ${parallel_jobs} -O ./

    # Switch environment for dbCAN analysis
    micromamba activate run_dbcan4

    # Run dbCAN on each part in parallel
    ls all_protein_families.part*.faa | /opt/conda/envs/COMPARATIVE_ANNOTATION/bin/parallel -j ${parallel_jobs} run_dbcan {} protein --dia_cpu 4 --hmm_cpu 4 --out_dir {}_dir --db_dir ${params.db_baseDir}/dbCAN \
    --hmm_eval  ${params.CAZyme_hmm_eval} --hmm_cov ${params.CAZyme_hmm_cov}

    # Merge the results
    mkdir db_can_out

    for file in overview.txt dbcan-sub.hmm.out diamond.out hmmer.out; do
        # Get the header from the first file
        head -n 1 \$(ls all_protein_families.part*_dir/\$file | head -n 1) > db_can_out/\$file
        # Append all data rows (excluding headers) from all files
        for f in all_protein_families.part*_dir/\$file; do
            tail -n +2 \$f >> db_can_out/\$file
        done
    done
    
    cat all_protein_families.part*_dir/uniInput > db_can_out/uniInput
    # cat all_protein_families.part*_dir/overview.txt > db_can_out/overview.txt
    # cat all_protein_families.part*_dir/dbcan-sub.hmm.out > db_can_out/dbcan-sub.hmm.out
    # cat all_protein_families.part*_dir/diamond.out > db_can_out/diamond.out
    # cat all_protein_families.part*_dir/hmmer.out > db_can_out/hmmer.out


    # Process the merged results
    cut -f1,3 db_can_out/overview.txt | grep -Pv "\\t-\$" | awk -F'\\t' '{split(\$2,a,"+"); for (i in a) {print \$1 "\\t" a[i];}}' | \
    awk -F"\\t" '{split(\$2,a,"_"); \$2=a[1]; print \$1 "\\t" \$2}' | sed -e "s/(.*//" > dbcan_family.tsv

    python ${params.scripts_baseDir}/dbcan_result_parse.py

    # Clean up intermediate files
    rm -rf all_protein_families.part*.faa all_protein_families.part*_dir  
    """
}



process run_dbCAN_visualization {
    publishDir "${params.visualization_results}/dbCAN", mode: 'copy'
    
    input:
    path(dbcan_HMMER_count_gene_count_matrix)
    path(dbcan_HMMER_count_gene_PA_matrix)
    path(metadata_file)

    output:
    path("heatmap_dbCAN_gene_PA_metadata${params.metacol}th_col.pdf")
    path("dbCAN_interactive_gene_PA_metadata${params.metacol}th_col")
    path("heatmap_dbCAN_gene_count_metadata${params.metacol}th_col.pdf")
    path("dbCAN_interactive_gene_count_metadata${params.metacol}th_col")

    when: params.run_visualization

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate R432_environment
    Rscript ${params.scripts_baseDir}/dbcan_visualization.R -i ${dbcan_HMMER_count_gene_PA_matrix} \
    -m ${metadata_file} -mc ${params.metacol} -out_html dbCAN_interactive_gene_PA_metadata${params.metacol}th_col
    mv heatmap_dbCAN.pdf heatmap_dbCAN_gene_PA_metadata${params.metacol}th_col.pdf

    Rscript ${params.scripts_baseDir}/dbcan_visualization.R -i ${dbcan_HMMER_count_gene_count_matrix} \
    -m ${metadata_file} -mc ${params.metacol} -out_html dbCAN_interactive_gene_count_metadata${params.metacol}th_col
    mv heatmap_dbCAN.pdf heatmap_dbCAN_gene_count_metadata${params.metacol}th_col.pdf
    """
}

process run_dbCAN {
    publishDir "${params.outdir}/dbCAN" , mode : 'copy'
    //conda "$HOME/miniforge3/envs/COMPARATIVE_ANNOTATION"
    cpus params.cpus

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 

    output:
    path("dbcan_gene_count_matrix.csv")
    path("dbcan_hmmerfamily_count_matrix.csv")
    path("heatmap_dbCAN.pdf")
    path("dbCAN_interactive")

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate run_dbcan4
    export OMP_NUM_THREADS=${task.cpus}

    time run_dbcan ppanggolin_result/pan_genes/all_protein_families.faa  protein --dia_cpu ${task.cpus}  \
    --hmm_cpu ${task.cpus}  --out_dir db_can_out  --db_dir ${params.db_baseDir}/dbCAN

    cut -f1,3 db_can_out/overview.txt | grep -Pv "\\t-\$" |  awk -F'\\t' '{split(\$2,a,"+"); for (i in a) {print \$1 "\\t" a[i];}}' | \
    awk -F"\\t" '{split(\$2,a,"_"); \$2=a[1]; print \$1 "\\t" \$2}'  | sed -e "s/(.*//" > dbcan_family.tsv

    python ${params.scripts_baseDir}/dbcan_result_parse.py 
    micromamba activate R432_environment
    Rscript  ${params.scripts_baseDir}/dbcan_visualization.R -i dbcan_hmmerfamily_count_matrix.csv  \
    -m ${metadata_file}  -mc ${params.metacol} -out heatmap_dbCAN.pdf -out_html dbCAN_interactive
    """
}

process run_eggNOG {
    publishDir "${params.annotation_results}/eggNOG" , mode : 'copy'
    cpus params.cpus

    input:
    path(ppanggolin_dir)

    output:
    path("eggnog_mmseqs.emapper*")
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate eggnog-mapper2112
    emapper.py -m mmseqs --data_dir ${params.db_baseDir}/eggNOG5 --output eggnog_mmseqs  --cpu ${task.cpus} \
    -i  ppanggolin_result/pan_genes/all_protein_families.faa --itype proteins
    """
}

process run_scoary2 {
    publishDir "${params.visualization_results}/scoary2" , mode : 'copy'
    cpus 8

    input:
    path(ppanggolin_dir)
    //path("${params.metadata}")
    path metadata_file 
    output:
    path("scoary_out")
    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate scoary-2
    python ${params.scripts_baseDir}/make_scoary_binary_trait.py -i ${metadata_file} -m ${params.metacol} -s ${params.samplecol}
    # genome_traits.tsv is generated

    # Set MPLCONFIGDIR to a writable directory
    export MPLCONFIGDIR=\$PWD/mpl_config
    mkdir -p \$MPLCONFIGDIR

    # Set NUMBA_CACHE_DIR to a writable directory
    export NUMBA_CACHE_DIR=\$PWD/numba_cache
    mkdir -p \$NUMBA_CACHE_DIR
    #export CONFINT_DB=\$PWD/confint_cache
    #mkdir -p \$CONFINT_DB
    export CONFINT_DB=$PWD/CONFINT_DB

    scoary2 --genes ppanggolin_result/gene_count_matrix.tsv  \
    --gene-data-type "gene-count:\\t" --traits genome_traits.tsv --trait-data-type "binary:\\t" --n-permut 1000 \
    --n-cpus ${task.cpus} --outdir scoary_out
    """
}


process run_drep_dereplication {
    publishDir "${params.outdir}/drep", mode: 'copy'
    cpus params.cpus

    input:
    path(genome_dir)

    output:
    path "drep_output"
    path "dereplicated_genomes"
    path "subspecies_clusters.tsv"


    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    #micromamba activate drep350

      #-g ${genome_dir}/*.fa \

    # Run dRep dereplication
    dRep dereplicate drep_output \
      -g *.fa \
      -p ${task.cpus} \
      --ignoreGenomeQuality \
      --S_algorithm ${params.drep_algorithm} \
      -sa ${params.drep_ani} \
      -nc ${params.drep_cov}
      
    # Copy dereplicated genomes to a separate directory for easier access
    mkdir -p dereplicated_genomes
    cp drep_output/dereplicated_genomes/*.fa dereplicated_genomes/

    echo -e "genome\tsubspecies_cluster\tprimary_cluster" > subspecies_clusters.tsv
    awk -F, 'NR>1 {print \$1"\\t"\$2"\\t"\$6}' drep_output/data_tables/Cdb.csv >> subspecies_clusters.tsv

    """
}


process countGenomes {
    input:
    path(genomes_dir)
    
    output:
    stdout emit: count
    
    script:
    """
    find "${genomes_dir}" -name "*.fa" | wc -l
    """
}

process create_shiny_dashboard {
    publishDir "${params.outdir}/shiny_dashboard", mode: 'copy'

    input:
    path visualization_results_dir

    output:
    path "shiny_dashboard.html"

    script:
    """
    #!/usr/bin/env python3
    import os

    html_content = "<html><body>"
    html_content += "<h1>Shiny App Links</h1>"

    for root, dirs, files in os.walk("${visualization_results_dir}"):
        for file in files:
            if file == "htShiny.sh":
                app_path = os.path.join(root, file)
                app_name = os.path.basename(os.path.dirname(app_path))
                relative_path = os.path.relpath(app_path, "${params.outdir}")
                html_content += f'<p><a href="{relative_path}" target="_blank">{app_name} Shiny App</a></p>'

    html_content += "</body></html>"

    with open("shiny_dashboard.html", "w") as f:
        f.write(html_content)
    """
}
process run_multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus 1

    input:
    path visualization_results
    path shiny_dashboard

    output:
    path "*"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate RAWREAD_QC
    export MPLCONFIGDIR=~/.config/matplotlib
    mkdir -p ~/.config/matplotlib


    multiqc ${visualization_results}  -o ${params.outdir}/multiqc
    cp ${shiny_dashboard} ${params.outdir}/multiqc/
    echo '<p>See <a href="shiny_dashboard.html">Shiny Apps Dashboard</a></p>' >> ${params.outdir}/multiqc/multiqc_report.html
    """
}

process make_sequence_db {
    publishDir "${params.outdir}/sequence_db", mode: 'copy'
    cpus params.cpus

    input:
    path metadata_file
    path prokka_dir
    path ppanggolin_dir

    output:
    path "sequences.h5"
    path "logs"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    #micromamba activate COMPARATIVE_ANNOTATION

    mkdir -p logs
    mkdir -p prokka
    mkdir -p ppanggolin_need 
        find -L ./ -name "*.faa" | grep -v ppanggolin_result | xargs -I {} cp {} prokka/
    find -L ./ -name "*.ffn" | grep -v ppanggolin_result | xargs -I {} cp {} prokka/
    
    cp ${ppanggolin_dir}/gene_families.tsv ppanggolin_need/
    cp ${ppanggolin_dir}/genome_id_linking.tsv ppanggolin_need/

    
    Rscript /scratch/tools/microbiome_analysis/comparative_annotation/make_sequence_db.R \
      --input ${metadata_file} \
      --processors ${task.cpus} \
      --output ./ \
      --data ./
    
 #rm -rf prokka
 #rm -rf ppanggolin_need
    """
}


workflow {
    take:
    annotation_results

    main:
    if (runPpanggolin) {
        ppanggolin_result = run_ppanggolin(prokka_collect)
    } else {
        ppanggolin_result = Channel.fromPath("${params.outdir}/ppanggolin_result")
    }
    if (runKofamscan) {
        kofamscan_result = run_kofamscan_annotation(ppanggolin_result)
    } else {
        kofamscan_result = Channel.fromPath("${params.annotation_results}/kofamscan/ko_matrix.csv")
    }
    if (runVFDB) {
        vfdb_result = run_VFDB_annotation(ppanggolin_result)
        vfdb_gene_pa = vfdb_result.vfdb_gene_pa
        vfdb_gene_count = vfdb_result.vfdb_gene_count
    } else {
        vfdb_gene_pa = Channel.fromPath("${params.annotation_results}/VFDB/gene_PA_VFDB_added.csv")
        vfdb_gene_count = Channel.fromPath("${params.annotation_results}/VFDB/gene_count_VFDB_added.csv")
    }

    if (runDbCAN) {
        dbcan_result = run_dbCAN_annotation(ppanggolin_result)
        dbcan_count_matrix = dbcan_result.HMMER_count_gene_count_matrix
        dbcan_pa_matrix = dbcan_result.HMMER_count_gene_PA_matrix
    } else {
        dbcan_count_matrix = Channel.fromPath("${params.annotation_results}/dbCAN/dbcan_HMMER_count_gene_count_matrix.csv")
        dbcan_pa_matrix = Channel.fromPath("${params.annotation_results}/dbCAN/dbcan_HMMER_count_gene_PA_matrix.csv")
    }
    if (runSkani) {
        all_genomes = selected_genomes.flatMap { dir -> file("${dir}/*.fa") }.collect()
        skani_result = run_skani_annotation(all_genomes)
        skani_fullmatrix = skani_result[0]
        skani_dist_result = skani_result[1]
    } else {
        skani_fullmatrix = Channel.fromPath("${params.annotation_results}/ani/skani_fullmatrix")
        skani_dist_result = Channel.fromPath("${params.annotation_results}/ani/skani_ANI_dist.tsv")
    }

    if (ca_skip_EggNOG) {
        eggnog_result = run_eggNOG(ppanggolin_result)
    }
    if (runMetadataSummary) {
        create_metadata_summary(metadata_ch)
    }
    if (run_sequence_db) {
        sequence_db = make_sequence_db(metadata_ch, prokka_collect, ppanggolin_result)
    } else {
        sequence_db = Channel.fromPath("${params.outdir}/sequence_db/sequences.h5")
    }

    if (runGenePACluster) {
        //if (fileCount >= 5) {
      //      log.info "Running GenePA Cluster with ${fileCount} genomes"
            run_genePA_cluster(ppanggolin_result, metadata_ch)
        //} else {
        //    log.info "Genome count (${fileCount}) is less than 5. GenePA Cluster will not be run."
        //}
    } else {
        log.info "GenePA Cluster already exists."
    }



         if (runDrep) {
             all_genomes = selected_genomes.flatMap { dir -> file("${dir}/*.fa") }.collect()
             drep_result = run_drep_dereplication(all_genomes)
         }


            metadata_ch = Channel.fromPath("${params.metadata}")


    if (params.run_visualization) {

            kofamscan_vis_result = run_kofamscan_visualization(
        kofamscan_result[0],
        metadata_ch          
    )
            //vfdb_vis_result = run_VFDB_visualization(vfdb_result.gene_pa, metadata_ch)
            vfdb_vis_result = run_VFDB_visualization(vfdb_gene_pa, vfdb_gene_count, metadata_ch)
            //card_vis_result = run_rgi_CARD_visualization(card_result.gene_pa, metadata_ch)
            card_vis_result = run_rgi_CARD_visualization(card_gene_pa, card_gene_count, metadata_ch)     
            dbcan_vis_result = run_dbCAN_visualization(dbcan_count_matrix, dbcan_pa_matrix, metadata_ch)
            //skani_vis_result = run_skani_visualization(skani_result[0], metadata_ch)
            skani_vis_result = run_skani_visualization(skani_fullmatrix, metadata_ch)
            //defensefinder_vis_result = run_defensefinder_visualization(defensefinder_gene_pa, defensefinder_gene_count, metadata_ch)

            scoary2_result = run_scoary2(ppanggolin_result, metadata_ch)

            vis_results = Channel.empty()
                .mix(kofamscan_vis_result)
                .mix(vfdb_vis_result)
                .mix(card_vis_result)
                .mix(dbcan_vis_result)
                .mix(skani_vis_result)
                .mix(scoary2_result)
                //.mix(defensefinder_vis_result) 
                .collect()
    } else {
        log.info "metacol is not specified. Visualization step will be skipped."
    }
            // Shiny dashboard creation and MultiQC execution after all visualizations are done
            //shiny_dashboard = create_shiny_dashboard(vis_results)
            //run_multiqc(vis_results, shiny_dashboard)


    }
