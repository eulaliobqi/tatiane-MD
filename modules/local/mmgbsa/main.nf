// MM-GBSA (trajetoria unica) via gmx_MMPBSA, sistema CHARMM36m+CGenFF
// (receptor 2I9T + daidzeina/UNL). Complementa (nao substitui) as analises
// GROMACS-nativas de ANALYSES/ANALYSES_RESIDUES — energia de ligacao com
// decomposicao por residuo, alvo direto dos residuos Arg30/Glu279 do docking.
//
// ── Licao do projeto irmao Milena-MD (modules/local/mmgbsa_robust/main.nf,
//    "abandonado" apos 3 tentativas) ────────────────────────────────────────
// A causa raiz mais provavel (nunca confirmada la, confirmada aqui por
// inspecao cruzada dos dois repos) NAO era "mamba run nao repassa argumentos
// posicionais" em geral -- e sim INDIRECAO DUPLA DE AMBIENTE CONDA:
//   1. O nextflow.config do Milena-MD NUNCA tinha um `withName: 'MMGBSA_ROBUST'
//      { conda = mmgbsa-env }` -- a task inteira rodava dentro do ambiente
//      DEFAULT (md-gromacs, GROMACS CUDA + AmberTools conflitantes), e o
//      modulo tentava "alcancar" mmgbsa-env de dentro do script via
//      `mamba run -n mmgbsa-env bash run_mmgbsa.sh "$1" "$2" "$3"`.
//   2. Rodar `mamba run -n <env>` a partir de dentro de uma task cujo shell
//      JA foi ativado por outro mecanismo de ativacao conda (o
//      `process.conda` do Nextflow tipicamente ativa o env via um script
//      wrapper `conda activate` antes de executar .command.sh) e uma
//      ativacao ANINHADA -- exatamente o cenario onde bugs de propagacao de
//      argv/env de conda/mamba run sao mais reportados (reconstroem o
//      ambiente do zero, inclusive PATH e as vezes $@ do subprocesso
//      intermediario). Isso explica por que MESMO a 3a tentativa (argumentos
//      posicionais, que deveriam atravessar o exec) ainda chegou com
//      -cs/-ct/-ci vazios: a indirecao dupla (ativacao do Nextflow + mamba
//      run manual por cima) e o problema estrutural, nao so a forma de
//      passar o argumento.
//   3. Faltava -cp (topologia do complexo) na chamada — bug secundario
//      independente, mascarado pelo primeiro erro.
//
// Fix aplicado aqui: este projeto (ver nextflow.config, profiles.conda.process)
// JA registra `withName: 'MMGBSA' { conda = mmgbsa-env }` — o Nextflow ativa
// mmgbsa-env nativamente para ESTA task, exatamente como faz para
// md-gromacs nos outros modulos (${params.gmx_cmd} chamado direto, nunca via
// `mamba run` dentro do script). Portanto `gmx_MMPBSA` e chamado DIRETO,
// sem nenhuma camada de ativacao/indirecao manual dentro do script — elimina
// a superficie do bug por completo. Todos os caminhos de arquivo
// (${md_tpr}, ${md_fit_xtc}, ${lig_ndx}, topol.top) sao interpolados pelo
// Groovy como texto literal ANTES do bash rodar; nao ha nenhuma "passagem de
// argumento entre processos/scripts" para dar errado.
process MMGBSA {
    tag "${meta.id}"
    label 'process_mmgbsa'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx),
                     path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*'),
                     path(prms, stageAs: 'prm_in/*')

    output:
    tuple val(meta), path("FINAL_RESULTS_MMGBSA.dat"), path("mmgbsa_results.csv"),
                     path("decomp_results.csv"), emit: results, optional: true
    tuple val(meta), path("mmgbsa.log"), path("mmgbsa_validation.txt"), emit: logs, optional: true

    script:
    def sys_name   = meta.id
    def saltcon    = params.nacl_conc        ?: 0.15
    // Frames: nao assumimos o dt de escrita da trajetoria (nao temos como
    // consultar isso aqui sem GROMACS confiavel no env) -- endframe bem
    // acima de qualquer trajetoria real e seguro (gmx_MMPBSA/cpptraj para no
    // ultimo frame disponivel, documentado). O "interval" e o unico knob que
    // realmente controla o tempo de execucao; exposto via params para poder
    // ajustar apos ver o numero real de frames no primeiro log (ver resumo).
    def startframe = params.mmgbsa_startframe ?: 1
    def endframe   = params.mmgbsa_endframe   ?: 9999999
    def interval   = params.mmgbsa_interval   ?: 10
    """
    echo "=== MMGBSA: ${meta.id} (gmx_MMPBSA, CHARMM36m+CGenFF, trajetoria unica) ===" >&2

    {
        echo "=== Validacao pre-MMGBSA: ${meta.id} ==="
        echo "Data: \$(date)"
        echo ""
    } > mmgbsa_validation.txt
    # mmgbsa.log criado vazio desde ja — se uma validacao antecipada abortar
    # com exit 0 antes do gmx_MMPBSA rodar, a tupla "logs" (mmgbsa.log +
    # mmgbsa_validation.txt) ainda assim e publicada (outputs "optional:
    # true" exigem TODOS os arquivos da tupla presentes p/ emitir; sem isso,
    # uma falha de validacao cedo faria a tupla inteira sumir).
    touch mmgbsa.log

    # ── 1. Materializa topologia + includes exatamente como TOPOLOGY/
    #        BOX_SOLVATE_IONS fazem (mesma convencao do resto do pipeline —
    #        NAO reinventar layout de #include) ──────────────────────────────
    cp input.top topol.top
    cp itp_in/*.itp . 2>/dev/null || true
    cp prm_in/*.prm . 2>/dev/null || true
    cp -r ${projectDir}/ff/charmm36-mar2019.ff .

    # Sanity check (nao fatal): confere que e a topologia JA solvatada
    # (BOX_SOLVATE_IONS), nao a pre-solvatacao (TOPOLOGY) — regra
    # nao-negociavel do laboratorio ("numero de coordenadas != topologia").
    # SOL/TIP3 = agua, NA/SOD = cation — nomenclatura CHARMM36 varia por porte.
    if ! grep -qiE "^\\s*(SOL|TIP3|HOH)\\s" topol.top; then
        echo "AVISO: 'SOL/TIP3/HOH' nao encontrado em [ molecules ] de topol.top — confirme que a topologia usada e a de BOX_SOLVATE_IONS (pos-solvatacao), nao a de TOPOLOGY (pre-solvatacao)" | tee -a mmgbsa_validation.txt
    fi

    # ── 2. Confirma grupos Receptor/Ligante em lig.ndx ───────────────────────
    if ! grep -q "\\[ *Receptor *\\]" ${lig_ndx}; then
        echo "ERRO FATAL: grupo 'Receptor' nao encontrado em ${lig_ndx}" | tee -a mmgbsa_validation.txt
        echo "No results — grupo Receptor ausente em lig.ndx" > FINAL_RESULTS_MMGBSA.dat
        echo "frame,TOTAL" > mmgbsa_results.csv
        echo "resid,resname,total" > decomp_results.csv
        exit 0
    fi
    if ! grep -q "\\[ *Ligante *\\]" ${lig_ndx}; then
        echo "ERRO FATAL: grupo 'Ligante' nao encontrado em ${lig_ndx}" | tee -a mmgbsa_validation.txt
        echo "No results — grupo Ligante ausente em lig.ndx" > FINAL_RESULTS_MMGBSA.dat
        echo "frame,TOTAL" > mmgbsa_results.csv
        echo "resid,resname,total" > decomp_results.csv
        exit 0
    fi
    echo "OK: grupos 'Receptor' e 'Ligante' confirmados em lig.ndx" >> mmgbsa_validation.txt

    # ── 3. Indice NUMERICO dos grupos (obrigatorio p/ -cg) ───────────────────
    # gmx_MMPBSA usa "-cg <indice_receptor> <indice_ligante>" com NUMEROS
    # (convencao GROMACS: blocos do .ndx numerados 0..N na ordem em que
    # aparecem no arquivo) -- NAO aceita nomes de grupo como "Receptor
    # Ligante" (esse era outro bug latente do modulo antigo do Milena-MD,
    # mascarado pelo erro de -cs/-ct/-ci vazio antes de chegar a este ponto;
    # confirmado na documentacao oficial: "-cg <Receptor group> <Ligand
    # group>", ie. -cg 1 13). Calculamos os indices aqui, na ordem real de
    # ${lig_ndx}, em vez de supor um numero fixo.
    IDX=0
    REC_IDX=""
    LIG_IDX=""
    while IFS= read -r line; do
        case "\$line" in
            \\[*\\])
                name=\$(echo "\$line" | sed -E 's/^\\[ *//; s/ *\\]\$//')
                [ "\$name" = "Receptor" ] && REC_IDX=\$IDX
                [ "\$name" = "Ligante" ]  && LIG_IDX=\$IDX
                IDX=\$((IDX + 1))
                ;;
        esac
    done < ${lig_ndx}

    echo "Indices calculados: Receptor=\${REC_IDX} Ligante=\${LIG_IDX} (grupos totais no .ndx: \${IDX})" | tee -a mmgbsa_validation.txt

    if [ -z "\${REC_IDX}" ] || [ -z "\${LIG_IDX}" ]; then
        echo "ERRO FATAL: nao foi possivel determinar indice numerico de Receptor/Ligante" | tee -a mmgbsa_validation.txt
        echo "No results — indice de grupo Receptor/Ligante nao determinado" > FINAL_RESULTS_MMGBSA.dat
        echo "frame,TOTAL" > mmgbsa_results.csv
        echo "resid,resname,total" > decomp_results.csv
        exit 0
    fi

    # ── 4. Pre-flight: gmx_MMPBSA precisa de um binario GROMACS (gmx/gmx_mpi)
    #        no PATH em tempo de execucao (usa editconf/trjconv internamente
    #        p/ preparar arquivos AMBER a partir do .tpr/.top) — mmgbsa-env
    #        foi desenhado so com AmberTools p/ nao colidir com o build CUDA
    #        de md-gromacs; se gromacs (CPU-only) nao foi instalado nesse env,
    #        falha aqui com mensagem clara em vez de um erro interno obscuro
    #        do gmx_MMPBSA la na frente. ────────────────────────────────────
    GMX_BIN=""
    for cand in gmx gmx_mpi; do
        command -v "\${cand}" >/dev/null 2>&1 && GMX_BIN="\${cand}" && break
    done
    if [ -z "\${GMX_BIN}" ]; then
        echo "AVISO: nenhum binario GROMACS (gmx/gmx_mpi) encontrado no PATH de mmgbsa-env — gmx_MMPBSA pode falhar internamente (editconf/trjconv). Ver resumo/comentarios deste modulo — pode ser necessario 'mamba install -n mmgbsa-env -c conda-forge gromacs' (build CPU-only, sem CUDA, p/ nao colidir com md-gromacs)." | tee -a mmgbsa_validation.txt
    else
        echo "OK: binario GROMACS encontrado em mmgbsa-env: \${GMX_BIN}" >> mmgbsa_validation.txt
    fi

    # ── 5. mmgbsa.in ──────────────────────────────────────────────────────
    # radiopt=0: CHARMM/CGenFF exige usar os raios ja definidos na propria
    # topologia (nao o mbondi2 default de sistemas AMBER) — documentado
    # oficialmente como obrigatorio p/ conversao GROMACS-CHARMM via ParmEd.
    # print_res="within 6": limita a decomposicao aos residuos na interface
    # (6 A do ligante) -- sem isso, idecomp=2 decompoe TODO residuo do
    # receptor (~275 no dominio de ligacao a DNA do 2I9T), inviabilizando o
    # tempo de execucao sem necessidade (nosso interesse sao os residuos de
    # interface, nao o receptor inteiro).
    cat > mmgbsa.in << MEOF
&general
sys_name="${sys_name}",
startframe=${startframe},
endframe=${endframe},
interval=${interval},
verbose=2,
radiopt=0,
/
&gb
igb=2,
saltcon=${saltcon},
/
&decomp
idecomp=2,
dec_verbose=1,
print_res="within 6",
/
MEOF
    echo "mmgbsa.in gerado (startframe=${startframe} endframe=${endframe} interval=${interval})" >> mmgbsa_validation.txt

    # ── 6. Executa gmx_MMPBSA — chamada DIRETA, sem wrapper/heredoc externo,
    #        sem mamba run (o Nextflow ja ativou mmgbsa-env nativamente via
    #        process.conda para esta task — ver nextflow.config). MPI
    #        paralelo (1 rank por frame) so se mpi4py estiver disponivel no
    #        env; senao roda serial e avisa (nao assumir MPI presente —
    #        environment-mmgbsa.yml historico so lista ambertools+gmx_MMPBSA,
    #        nao mpi4py explicitamente). ──────────────────────────────────
    MPI_PREFIX=""
    if command -v mpirun >/dev/null 2>&1 && python3 -c "import mpi4py" >/dev/null 2>&1; then
        MPI_PREFIX="mpirun -np ${task.cpus}"
        echo "MPI disponivel — rodando com \${MPI_PREFIX}" >> mmgbsa_validation.txt
    else
        echo "mpi4py/mpirun nao encontrado em mmgbsa-env — rodando gmx_MMPBSA serial (mais lento; considere 'mamba install -n mmgbsa-env -c conda-forge mpi4py openmpi')" | tee -a mmgbsa_validation.txt
    fi

    echo "[MMGBSA] Iniciando gmx_MMPBSA (pode demorar — ate 72h de budget)..." >&2

    \${MPI_PREFIX} gmx_MMPBSA -O \\
        -i   mmgbsa.in \\
        -cs  ${md_tpr} \\
        -ct  ${md_fit_xtc} \\
        -ci  ${lig_ndx} \\
        -cg  \${REC_IDX} \${LIG_IDX} \\
        -cp  topol.top \\
        -o   FINAL_RESULTS_MMGBSA.dat \\
        -eo  mmgbsa_results.csv \\
        -deo decomp_results.csv \\
        -nogui \\
        2>&1 | tee mmgbsa.log

    # ── 7. Verifica saidas e cria fallbacks graciosos ────────────────────────
    # errorStrategy 'ignore' + outputs optional: true => se a task terminar
    # com exit != 0, NENHUM output e publicado (nem os placeholders escritos
    # no work dir) -- por isso sempre terminamos com exit 0, mesmo em falha,
    # exatamente como o modulo antigo do Milena-MD ja fazia (e documentou
    # como proposital). REPORT deste projeto nao depende do canal de saida
    # do MMGBSA (le o resultado do disco), entao o placeholder aqui e so
    # para diagnostico humano, nunca para desbloquear o resto do pipeline.
    if [ -s FINAL_RESULTS_MMGBSA.dat ] && ! grep -q "^No results" FINAL_RESULTS_MMGBSA.dat; then
        echo "[OK] FINAL_RESULTS_MMGBSA.dat gerado" | tee -a mmgbsa_validation.txt
        [ -s decomp_results.csv ] || {
            echo "resid,resname,total" > decomp_results.csv
            echo "AVISO: decomp_results.csv nao gerado — arquivo vazio criado" >> mmgbsa_validation.txt
        }
        [ -s mmgbsa_results.csv ] || {
            echo "frame,TOTAL" > mmgbsa_results.csv
            echo "AVISO: mmgbsa_results.csv nao gerado — arquivo vazio criado" >> mmgbsa_validation.txt
        }
        echo "[MMGBSA] Concluido com sucesso para ${meta.id}" >&2
    else
        echo "ERRO: FINAL_RESULTS_MMGBSA.dat nao gerado (ou vazio)" | tee -a mmgbsa_validation.txt
        echo "--- Ultimas 40 linhas do log ---" >> mmgbsa_validation.txt
        tail -40 mmgbsa.log >> mmgbsa_validation.txt 2>/dev/null || true
        echo "No results — gmx_MMPBSA failed (ver mmgbsa.log / mmgbsa_validation.txt)" > FINAL_RESULTS_MMGBSA.dat
        echo "frame,TOTAL" > mmgbsa_results.csv
        echo "resid,resname,total" > decomp_results.csv
        echo "[MMGBSA] FALHOU para ${meta.id} — ver mmgbsa_validation.txt" >&2
    fi
    """
}
