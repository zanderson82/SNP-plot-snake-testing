rule make_SNP_plot:
    input:
        filtered_files = f"{config['outputpath']}/filtered-annotated-VCFs/file_list.txt"
    output:
        plots = f"{config['outputpath']}/SNP-Plots/{config['description']}-SNP-PLOT.svg"
    params:
        demographic_data = config["demographicFile"]
    conda:
        "Rtools-1.1"
    script:
        "/n/zanderson/SNP-plot-snake-testing/workflow/scripts/SNP-plot-1.R"