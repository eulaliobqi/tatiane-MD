// Variante de ANALYSES para molecula pequena (UNL/daidzeina), mirror de
// MD-gromacs/modules/local/analyses_ben — detecta cadeia B via HETATM (nao
// ATOM) para encontrar o ligante. ANALYSES_RESIDUES nao precisa alteracao —
// recebe lig.ndx com "Ligante"/"Receptor" ja definidos.
process ANALYSES {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(md_tpr), path(md_fit_xtc)

    output:
    tuple val(meta), path("*.xvg"), path("lig.ndx"), emit: xvg

    script:
    """
    # UNL e HETATM na cadeia B — detectar residuo correto
    LIG_FIRST=\$(awk '(/^ATOM/ || /^HETATM/) && substr(\$0,22,1)=="B" {print substr(\$0,23,4)+0; exit}' ${complexo_pdb})
    LIG_LAST=\$(awk '(/^ATOM/ || /^HETATM/) && substr(\$0,22,1)=="B" {r=substr(\$0,23,4)+0} END{print r}' ${complexo_pdb})

    if [ -z "\${LIG_FIRST}" ] || [ -z "\${LIG_LAST}" ]; then
        echo "ERRO: cadeia B (UNL) nao detectada em ${complexo_pdb}"; exit 1
    fi
    echo "[ANALISES] Ligante UNL: residuo \${LIG_FIRST}-\${LIG_LAST}" >&2

    N_DEFAULT=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -o _default.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    LIG_IDX=\${N_DEFAULT}
    REC_IDX=\$((N_DEFAULT + 1))
    rm -f _default.ndx

    ${params.gmx_cmd} make_ndx -f ${md_tpr} -o lig.ndx << EOF
r \${LIG_FIRST}-\${LIG_LAST}
name \${LIG_IDX} Ligante
1 & ! \${LIG_IDX}
name \${REC_IDX} Receptor
q
EOF

    # RMSD backbone do receptor
    printf 'Backbone\\nBackbone\\n' | ${params.gmx_cmd} rms \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o rmsd_backbone.xvg -tu ns

    # RMSD do ligante
    printf 'Ligante\\nLigante\\n' | ${params.gmx_cmd} rms \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o rmsd_ligante.xvg -tu ns

    # RMSF por residuo (backbone receptor)
    printf 'Backbone\\n' | ${params.gmx_cmd} rmsf \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o rmsf_residuos.xvg -res -fit

    # Raio de giro (proteina)
    printf 'Protein\\n' | ${params.gmx_cmd} gyrate \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o gyrate.xvg -tu ns

    # Contatos receptor-ligante < 0.4 nm
    printf 'Receptor\\nLigante\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -od mindist_receptor_ligante.xvg -on numcont_receptor_ligante.xvg \\
        -d 0.4 -tu ns

    # Pontes de hidrogenio receptor-ligante — CGenFF usa nomes de atomo CHARMM
    # padrao (N/O+H), gmx hbond costuma reconhecer normalmente; mantem
    # fallback defensivo (mesmo padrao usado no pipeline BEN/GAFF2, onde isso
    # de fato falhou em producao naquele caso)
    printf 'Receptor\\nLigante\\n' | ${params.gmx_cmd} hbond \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -num hbond.xvg -tu ns 2>&1 | tee hbond.log || \\
        printf '# gmx hbond falhou — ver hbond.log\\n@ title "Number of Hydrogen Bonds"\\n@ xaxis label "Time (ns)"\\n@ yaxis label "Number"\\n0.000 0\\n' > hbond.xvg

    # SASA do receptor
    ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx \\
        -surface 'Protein' -output 'Protein' \\
        -o sasa_receptor.xvg -tu ns

    # SASA do ligante (valores baixos = UNL enterrado na interface)
    ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx \\
        -surface 'Ligante' -output 'Ligante' \\
        -o sasa_ligante.xvg -tu ns

    echo "[OK] Analises concluidas para ${meta.id}" >&2
    """
}
