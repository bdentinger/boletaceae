# Build the static site regardless of whether a release exists yet.
# build_site.py already handles empty/partial states.

rule site_build:
    output:
        directory("site")
    shell:
        "{config[python_bin]} workflow/scripts/build_site.py && touch site/.nojekyll"

# Keep rule all minimal here if you have another rule all in the root Snakefile.
