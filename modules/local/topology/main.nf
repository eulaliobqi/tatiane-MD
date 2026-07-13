process TOPOLOGY {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/topo" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb),
                     path(unl_itp,        stageAs: 'lig_in/unl.itp'),
                     path(unl_prm,        stageAs: 'lig_in/unl.prm'),
                     path(posre_unl_itp,  stageAs: 'lig_in/posre_UNL.itp'),
                     path(unl_ini_gro,    stageAs: 'lig_in/unl_ini.gro')

    output:
    tuple val(meta), path("complexo.gro"), path("topol.top"), path("*.itp"), path("*.prm"),
                     emit: topology

    script:
    """
    echo "=== TOPOLOGY: ${meta.id} (pdb2gmx charmm36-feb2026_cgenff-5.0 + merge CGenFF) ===" >&2

    # gmx pdb2gmx localiza NAME.ff por nome apenas se a pasta existir no cwd
    cp -r ${projectDir}/ff/charmm36-feb2026_cgenff-5.0.ff .

    awk '/^ATOM/ && substr(\$0,22,1)=="A" {print}' ${complexo_pdb} > receptor.pdb
    echo "TER" >> receptor.pdb
    echo "END" >> receptor.pdb

    NATOM_REC=\$(grep -c "^ATOM" receptor.pdb || echo 0)
    echo "  Receptor: \${NATOM_REC} atomos ATOM (cadeia A)" >&2

    # Sem -ter: termos padrao (NH3+/COO-) sem prompt interativo — os indices
    # numericos do prompt -ter diferem entre AMBER e CHARMM, nao reusar o
    # "printf '0\\n0\\n'" do pipeline BEN/AMBER aqui.
    # -ignh e redundante (receptor ja sem H apos pdb2pqr) mas mantido por seguranca.
    ${params.gmx_cmd} pdb2gmx \\
        -f receptor.pdb \\
        -o receptor.gro \\
        -p receptor.top \\
        -i posre.itp \\
        -ff charmm36-feb2026_cgenff-5.0 \\
        -water tip3p \\
        -ignh \\
        2>&1 | tee pdb2gmx.log

    if [ ! -s receptor.gro ]; then
        echo "ERRO: pdb2gmx falhou — ver pdb2gmx.log"; exit 1
    fi

    python3 ${projectDir}/bin/merge_small_molecule_topology.py \\
        --protein-gro receptor.gro \\
        --ligand-gro  ${unl_ini_gro} \\
        --protein-top receptor.top \\
        --ligand-itp  ${unl_itp} \\
        --ligand-prm  ${unl_prm} \\
        --ligand-mol  UNL \\
        --out-gro     complexo.gro \\
        --out-top     topol.top

    if [ ! -s complexo.gro ]; then
        echo "ERRO: merge_small_molecule_topology.py falhou"; exit 1
    fi

    # Copia os arquivos do ligante do subdiretorio de staging (lig_in/) para o
    # cwd — o Nextflow EXCLUI arquivos de input do casamento de glob de output
    # por padrao ("input files are not included in the default matching set"),
    # entao path("*.itp")/path("*.prm") nao enxergariam unl.itp/unl.prm/
    # posre_UNL.itp se eles ficassem só no local onde foram staged como input.
    # Pego rodando de verdade no servidor (Missing output file(s) `*.prm`).
    cp lig_in/unl.itp .
    cp lig_in/unl.prm .
    cp lig_in/posre_UNL.itp .

    NTOTAL=\$(awk 'NR==2{print \$1}' complexo.gro)
    echo "[OK] TOPOLOGY concluido: complexo.gro com \${NTOTAL} atomos" >&2
    """
}
