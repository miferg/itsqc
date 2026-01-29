import subprocess
import glob
import os
import pathlib

querydir = Path(config["querydir"])
outdir = Path(config["outdir"])
LOCBASE = [x.split('.fna')[0] for x in os.listdir(querydir) if x.endswith('.fna')]

rule all:
    input:
        str(outdir) + "/merged/itsqc.merged.fna",
        str(outdir) + "/clean/itsqc.clean.fna",
        str(outdir) + "/clean/itsqc.clean.o.fna",
        str(outdir) + "/clean/itsqc.derep.fna",
        str(outdir) + "/itsx_out/itsx_out.positions.txt",
        str(outdir) + "/itsqc.wits.tsv",
        str(outdir) + "/itsqc.wits.fna",
        str(outdir) + "/annot/sh_general_release_dynamic_s_19.02.2025.fasta"

rule merge_fna_files:
    input:
        expand(str(querydir) + "/{lbase}.fna", lbase=LOCBASE)
    output:
        str(outdir) + "/merged/itsqc.merged.fna"
    params:
        str(querydir) +"/*.fna"
    shell:
        """
        cat {params} > {output}
        """

rule clean_sequences:
    conda:
        "snakes/itsqc_pylibs.yaml"
    input:
        str(outdir) + "/merged/itsqc.merged.fna"
    output:
        str(outdir) + "/clean/itsqc.clean.fna"
    params:
        str(outdir) + "/clean/itsqc"
    shell:
        """
        python snakes/trim_and_filter.py {input} {params}
        """

rule get_annot_db:
   # download reference database
    output:
        temp(str(outdir) + "/annot/General_EUK_ITS_v2.0.zip")
    params:
        str(outdir) +'/annot'
    shell:
        """
        wget --directory-prefix={params} https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v2.0.zip
        """
    
rule unzip_annot_db:
   # decompress reference database
    input:
        str(outdir) + "/annot/General_EUK_ITS_v2.0.zip"
    output:
        str(outdir) + "/annot/General_EUK_ITS_v2.0.fasta"
    params:
        str(outdir) +'/annot'
    shell:
        """
        7zip x {input} -o{params}
        """

rule  orient_sequences:
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/clean/itsqc.clean.fna",
        str(outdir) + "/annot/General_EUK_ITS_v2.0.fasta"
    output:
        str(outdir) + "/clean/itsqc.clean.o.fna"
    shell:
        """
        vsearch --orient {input[0]} --db {input[1]} --fastaout {output}
        """

rule dereplicate:
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/clean/itsqc.clean.o.fna"
    output:
        str(outdir) + "/clean/itsqc.derep.fna",
        str(outdir) + "/clean/itsqc.derep.uc"
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --uc {output[1]}
        """

rule search_itsx:
    conda:
        "snakes/itsx.yaml"
    input:
        str(outdir) + "/clean/itsqc.derep.fna"
    output:
        str(outdir) + "/itsx_out/itsx_out.positions.txt"
    params:
        str(outdir) + "/itsx_out/itsx_out"
    shell:
        """
        ITSx -i {input} -o {params} --cpu 4 -t all
        """

rule filter_itsx_table:
    input:
        str(outdir) + "/itsx_out/itsx_out.positions.txt"
    output:
        str(outdir) + "/itsqc.wits.tsv" # with ITS
    shell:
        """
        awk -F"\t" '{{print $1"\t"$4"\t"$6}}' {input} | grep -v "Not found" > {output}
        """

rule  filter_by_its_content:
    conda:
        "snakes/itsqc_pylibs.yaml"
    input:
        str(outdir) + "/clean/itsqc.derep.fna",
        str(outdir) + "/itsqc.wits.tsv"
    output:
        str(outdir) + "/itsqc.wits.fna"
    params:
        str(outdir) + "/itsqc"
    shell:
        """
        python snakes/filter_by_its_content.py {input[0]} {input[1]} {params}
        """

# END
