rule site_build:
    input:
        expand("data/releases/{v}/manifest.json", v=[config["version"]])
    output:
        directory("site")
    shell:
        "python3 workflow/scripts/build_site.py"
