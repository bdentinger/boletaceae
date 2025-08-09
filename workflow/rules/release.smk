rule site_build:
    input:
        expand("data/releases/{v}/manifest.json", v=[config["version"]])
    output:
        directory("site")
    shell:
        "{config[python_bin]} workflow/scripts/build_site.py"
