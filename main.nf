#!/usr/bin/env nextflow
// Pipeline de Dinamica Molecular — receptor 2I9T (dominio de ligacao a DNA de
// NF-kB, cadeia A) + ligante daidzeina (isoflavona, resname UNL, pose de
// docking AutoDock Vina). CHARMM36m + CGenFF 5.0 (ParamChem) — distinto do
// AMBER99SB-ILDN+GAFF2/ACPYPE usado nos demais pipelines deste laboratorio
// (MD-gromacs, Milena-MD). Ver README.md para as decisoes de parametros
// (pH 7.4, NaCl 0.15M, box 1.2nm — alvo humano, nao inseto).
nextflow.enable.dsl = 2

include { PREPARE_PH          } from './modules/local/prepare_ph/main.nf'
include { LIGAND_TOPOLOGY     } from './modules/local/ligand_topology/main.nf'
include { PREPARE_COMPLEX     } from './modules/local/prepare_complex/main.nf'
include { TOPOLOGY            } from './modules/local/topology/main.nf'
include { BOX_SOLVATE_IONS    } from './modules/local/box_solvate_ions/main.nf'
include { MINIMIZATION        } from './modules/local/minimization/main.nf'
include { NVT                 } from './modules/local/nvt/main.nf'
include { NPT                 } from './modules/local/npt/main.nf'
include { PRODUCTION          } from './modules/local/production/main.nf'
include { POSTPROCESS         } from './modules/local/postprocess/main.nf'
include { ANALYSES            } from './modules/local/analyses/main.nf'
include { ANALYSES_RESIDUES   } from './modules/local/analyses_residues/main.nf'
include { MMGBSA              } from './modules/local/mmgbsa/main.nf'
include { PLOT                } from './modules/local/plot/main.nf'
include { REPORT              } from './modules/local/report/main.nf'

workflow {
    if (!params.input) {
        error "Informe o samplesheet: --input samplesheet.csv"
    }

    // Canal do receptor: (meta, receptor.pdb)
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta     = [id: row.sample_id]
            def receptor = file(row.receptor, checkIfExists: true)
            tuple(meta, receptor)
        }
        .set { ch_receptor }

    // Canal do ligante: (meta, mol2, str) — independente do receptor, roda em paralelo
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample_id]
            def mol2 = file(row.ligand_mol2, checkIfExists: true)
            def str  = file(row.ligand_str,  checkIfExists: true)
            tuple(meta, mol2, str)
        }
        .set { ch_ligand }

    PREPARE_PH(ch_receptor)
    LIGAND_TOPOLOGY(ch_ligand)

    // Junta por meta.id (String) — NAO por meta (Map) inteiro: Maps com
    // chaves diferentes nunca sao iguais em Groovy, join() por Map inteiro
    // ja silenciosamente nao emitiu nada em producao neste laboratorio
    // (ver Milena-MD/main.nf, comentario sobre ANALYSES_TRIAD). Aqui os
    // metas sao identicos (so {id:...}) mas mantemos a convencao por
    // seguranca e consistencia entre todos os joins do pipeline.
    ch_complex_input = PREPARE_PH.out.protonated
        .map { meta, receptor_ph -> tuple(meta.id, meta, receptor_ph) }
        .join(
            LIGAND_TOPOLOGY.out.topology
                .map { meta, itp, prm, posre, ini_pdb, ini_gro -> tuple(meta.id, ini_pdb) },
            by: 0
        )
        .map { id, meta, receptor_ph, ini_pdb -> tuple(meta, receptor_ph, ini_pdb) }

    PREPARE_COMPLEX(ch_complex_input)

    ch_topology_input = PREPARE_COMPLEX.out.complexo
        .map { meta, complex_pdb -> tuple(meta.id, meta, complex_pdb) }
        .join(
            LIGAND_TOPOLOGY.out.topology
                .map { meta, itp, prm, posre, ini_pdb, ini_gro -> tuple(meta.id, itp, prm, posre, ini_gro) },
            by: 0
        )
        .map { id, meta, complex_pdb, itp, prm, posre, ini_gro ->
            tuple(meta, complex_pdb, itp, prm, posre, ini_gro)
        }

    TOPOLOGY(ch_topology_input)

    BOX_SOLVATE_IONS(TOPOLOGY.out.topology)
    MINIMIZATION(BOX_SOLVATE_IONS.out.system)
    NVT(MINIMIZATION.out.system)
    NPT(NVT.out.system)
    PRODUCTION(NPT.out.system)
    POSTPROCESS(PRODUCTION.out.traj)

    // ANALYSES precisa do complex.pdb (deteccao de cadeia B/HETATM) + trajetoria pos-processada
    ch_analyses_input = PREPARE_COMPLEX.out.complexo
        .map { meta, complex_pdb -> tuple(meta.id, meta, complex_pdb) }
        .join(
            POSTPROCESS.out.fit.map { meta, tpr, xtc -> tuple(meta.id, tpr, xtc) },
            by: 0
        )
        .map { id, meta, complex_pdb, tpr, xtc -> tuple(meta, complex_pdb, tpr, xtc) }

    ANALYSES(ch_analyses_input)

    // ANALYSES_RESIDUES recebe lig.ndx de ANALYSES — funciona sem alteracao
    ch_residues_input = POSTPROCESS.out.fit
        .map { meta, tpr, xtc -> tuple(meta.id, meta, tpr, xtc) }
        .join(
            ANALYSES.out.xvg.map { meta, xvgs, ndx -> tuple(meta.id, ndx) },
            by: 0
        )
        .map { id, meta, tpr, xtc, ndx -> tuple(meta, tpr, xtc, ndx) }

    ANALYSES_RESIDUES(ch_residues_input)

    // MMGBSA usa a topologia JA solvatada/com ions (BOX_SOLVATE_IONS, NAO
    // TOPOLOGY) — regra nao-negociavel deste laboratorio: o .top pre-
    // solvatacao nao bate em contagem de atomos com a trajetoria de producao
    // (ver bioinformatics.md). errorStrategy 'ignore' no modulo: falha aqui
    // NAO deve derrubar o resto do pipeline (MM-GBSA e suplementar).
    ch_mmgbsa_input = POSTPROCESS.out.fit
        .map { meta, tpr, xtc -> tuple(meta.id, meta, tpr, xtc) }
        .join(
            ANALYSES.out.xvg.map { meta, xvgs, ndx -> tuple(meta.id, ndx) },
            by: 0
        )
        .join(
            BOX_SOLVATE_IONS.out.system
                .map { meta, gro, top, itps, prms -> tuple(meta.id, top, itps, prms) },
            by: 0
        )
        .map { id, meta, tpr, xtc, ndx, top, itps, prms ->
            tuple(meta, tpr, xtc, ndx, top, itps, prms)
        }

    MMGBSA(ch_mmgbsa_input)

    // PLOT e REPORT NAO dependem do canal de saida do MMGBSA (errorStrategy
    // 'ignore' pode nunca emitir; um join contra esse canal ja travou o PLOT
    // inteiro em producao no projeto irmao Milena-MD, mesmo com
    // remainder:true — ver comentario em modules/local/report/main.nf).
    // REPORT le o resultado do MMGBSA como caminho de disco (publishDir),
    // com fallback gracioso se ainda nao existir/tiver falhado.
    ch_plot_input = ANALYSES.out.xvg
        .map { meta, xvgs, ndx -> tuple(meta.id, meta, xvgs) }
        .join(
            ANALYSES_RESIDUES.out.residues
                .map { meta, d1, d2, s1, s2 -> tuple(meta.id, [d1, d2, s1, s2]) },
            by: 0
        )
        .map { id, meta, xvgs, res_files -> tuple(meta, xvgs, res_files) }

    PLOT(ch_plot_input)
    REPORT(ch_plot_input)
}
