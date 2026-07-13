process TOPOLOGY {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/topo" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(unl_itp), path(unl_prm),
                     path(posre_unl_itp), path(unl_ini_gro)

    output:
    tuple val(meta), path("complexo.gro"), path("topol.top"), path("*.itp"), path("*.prm"),
                     emit: topology

    script:
    """
    echo "=== TOPOLOGY: ${meta.id} (pdb2gmx charmm36-mar2019 + merge CGenFF) ===" >&2

    # gmx pdb2gmx localiza NAME.ff por nome apenas se a pasta existir no cwd
    cp -r ${projectDir}/ff/charmm36-mar2019.ff .

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
        -ff charmm36-mar2019 \\
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

    # posre_unl_itp ja esta no cwd (staged pelo Nextflow como input declarado) —
    # nao precisa copiar; um `cp` pra si mesmo aqui gerava erro espurio
    # ("are the same file"), pego rodando de verdade no servidor
    NTOTAL=\$(awk 'NR==2{print \$1}' complexo.gro)
    echo "[OK] TOPOLOGY concluido: complexo.gro com \${NTOTAL} atomos" >&2
    """
}
