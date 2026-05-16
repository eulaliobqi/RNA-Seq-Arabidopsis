Diagnostique falhas no pipeline Nextflow usando o agente pipeline-monitor.

Leia os seguintes arquivos de log (use os que existirem):
- .nextflow.log (log principal do Nextflow)
- results/nextflow_trace.txt (trace de processos)
- logs/ (logs de execução do run_pipeline.sh)

Se o pipeline falhou, identifique:
1. Qual processo falhou (nome do processo Nextflow)
2. Qual o erro exato (mensagem de erro)
3. A causa provável (memória, arquivo ausente, ferramenta, R, etc.)
4. A correção recomendada (arquivo + linha + valor)
5. Se é possível retomar com --resume

Seja direto e forneça o diagnóstico em no máximo 20 linhas.
