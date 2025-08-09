# (Optional) placeholder for future align/trim rules.
# Not required for the static site to build, but kept for structure.
# You can add MAFFT/TRIMAL steps later.

# Example stub rule to make Snakemake happy if referenced:
rule noop:
    output: "data/qc/.touch"
    shell: "mkdir -p data/qc && touch {output}"
