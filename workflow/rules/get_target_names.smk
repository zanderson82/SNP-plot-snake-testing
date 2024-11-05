rule get_target_names:
    input: 
        sample_dir = config["sampledirectory"],
        target_file = config["targetfile"]
    output:
        input_vcf_list = temp("/n/zanderson/SNP-plot-snake-testing/config/input_vcf_list.txt")
    conda:
        "iPython-8.15.0"
    script:
        "/n/zanderson/SNP-plot-snake-testing/workflow/scripts/get_target_names.py"

