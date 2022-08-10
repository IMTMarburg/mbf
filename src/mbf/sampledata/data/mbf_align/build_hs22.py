import mbf_genomes
import subprocess


genome = mbf_genomes.EnsemblGenome('Homo_sapiens', 94)
allowed = set(['22'])
genes = genome.df_genes[genome.df_genes['chr'].isin(allowed)]
transcripts = genome.df_transcripts[genome.df_transcripts['chr'].isin(allowed)]
genes.to_msgpack("hs_22_genes.msgpack")
transcripts.to_msgpack("hs_22_transcripts.msgpack")

subprocess.check_call(['gzip', 'hs_22_genes.msgpack'])
subprocess.check_call(['gzip', 'hs_22_transcripts.msgpack'])
