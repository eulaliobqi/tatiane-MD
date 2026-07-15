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

// Analises extras (2026-07-14, ver memoria do projeto) — portadas de
// Milena-MD/MD-gromacs, adaptadas pra ligante molecula-pequena. Rodam sobre
// as duas fases da trajetoria (bound 0-60ns / relocated 65ns-fim, ver
// PHASE_SPLIT), nao a trajetoria inteira — a transicao real detectada em
// ~65-78ns tornaria uma media/analise unica sobre os 100ns enganosa.
include { PHASE_SPLIT         } from './modules/local/phase_split/main.nf'
include { CLUSTERING          } from './modules/local/clustering/main.nf'
include { CONTACT_MAP         } from './modules/local/contact_map/main.nf'
include { PROLIF_FINGERPRINT  } from './modules/local/prolif_fingerprint/main.nf'
include { FE_RERUN            } from './modules/local/fe_rerun/main.nf'
include { FE_INTERPRET        } from './modules/local/fe_interpret/main.nf'

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

    // ══════════════════════════════════════════════════════════════════════
    // Analises extras (2026-07-14): mapa de contatos, fingerprint ProLIF,
    // clustering e estimativa de energia (Interaction Entropy) — todas
    // tolerantes a falha (errorStrategy 'ignore'), rodando sobre as DUAS
    // fases da trajetoria (bound/relocated), nao sobre os 100ns inteiros.
    // ══════════════════════════════════════════════════════════════════════

    // PHASE_SPLIT reusa tpr/xtc de POSTPROCESS + lig.ndx de ANALYSES —
    // mesma entrada que MMGBSA ja usa.
    ch_phase_split_input = POSTPROCESS.out.fit
        .map { meta, tpr, xtc -> tuple(meta.id, meta, tpr, xtc) }
        .join(
            ANALYSES.out.xvg.map { meta, xvgs, ndx -> tuple(meta.id, ndx) },
            by: 0
        )
        .map { id, meta, tpr, xtc, ndx -> tuple(meta, tpr, xtc, ndx) }

    PHASE_SPLIT(ch_phase_split_input)

    // Funde as duas fases num unico canal, tag meta.phase p/ diferenciar
    // publishDir/tag de cada task downstream (mesmo meta.id, phase distinto
    // — por isso todo join daqui pra baixo usa meta.id, nunca o Map inteiro).
    ch_phase_all = PHASE_SPLIT.out.bound
        .map { meta, tpr, xtc, ndx -> tuple(meta + [phase: 'bound'], tpr, xtc, ndx) }
        .mix(
            PHASE_SPLIT.out.relocated
                .map { meta, tpr, xtc, ndx -> tuple(meta + [phase: 'relocated'], tpr, xtc, ndx) }
        )

    CLUSTERING(ch_phase_all)

    // CONTACT_MAP/PROLIF_FINGERPRINT precisam do complexo.pdb original (nao
    // do .tpr/.gro) — mesmo motivo do Milena-MD: deteccao de cadeia B usa
    // leitura direta do PDB. complexo.pdb NAO e' por fase — broadcast do
    // MESMO valor pras 2 linhas de fase.
    //
    // NAO usar .join() aqui: join() pareia 1:1 e CONSOME o item do lado
    // direito -- com 2 linhas do lado esquerdo (bound/relocated) com a
    // MESMA chave (meta.id) contra 1 unico item de complexo.pdb, apenas a
    // PRIMEIRA linha a chegar casa; a segunda fica sem par e e' descartada
    // silenciosamente (confirmado em producao 2026-07-14: CONTACT_MAP/
    // PROLIF_FINGERPRINT/FE_RERUN/FE_INTERPRET rodaram so "1 of 1" em vez
    // de "2 of 2", cada um pegando uma fase diferente por corrida de
    // conclusao das tasks upstream). .combine(by:0) faz o produto
    // cartesiano correto por chave (2 esquerda x 1 direita = 2 saidas).
    ch_interact_input = ch_phase_all
        .map { meta, tpr, xtc, ndx -> tuple(meta.id, meta, tpr, xtc) }
        .combine(
            PREPARE_COMPLEX.out.complexo.map { meta, complex_pdb -> tuple(meta.id, complex_pdb) },
            by: 0
        )
        .map { id, meta, tpr, xtc, complex_pdb -> tuple(meta, complex_pdb, tpr, xtc) }

    CONTACT_MAP(ch_interact_input)
    PROLIF_FINGERPRINT(ch_interact_input)

    // FE_RERUN precisa da topologia pos-solvatacao (BOX_SOLVATE_IONS, NAO
    // TOPOLOGY — mesma regra nao-negociavel do MMGBSA) + estrutura final da
    // producao (PRODUCTION.out.checkpoint) + subtrajetoria reduzida de
    // CLUSTERING. Nem PRODUCTION nem BOX_SOLVATE_IONS sao por fase —
    // mesmo motivo acima, .combine(by:0) em vez de .join(by:0).
    ch_fe_input = CLUSTERING.out.for_fe
        .map { meta, tpr, xtc, ndx -> tuple(meta.id, meta, xtc, ndx) }
        .combine(
            PRODUCTION.out.checkpoint.map { meta, gro, cpt, edr -> tuple(meta.id, gro) },
            by: 0
        )
        .combine(
            BOX_SOLVATE_IONS.out.system
                .map { meta, gro, top, itps, prms -> tuple(meta.id, top, itps) },
            by: 0
        )
        .map { id, meta, xtc, ndx, prod_gro, top, itps ->
            tuple(meta, prod_gro, top, itps, xtc, ndx)
        }

    FE_RERUN(ch_fe_input)
    FE_INTERPRET(FE_RERUN.out.energy)

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
