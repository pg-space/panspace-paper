# !/bin/bash

# DIRSAVE="data/bacteria-all"
# DIRSAVE_TEST="data/bacteria-test"
DOWNLOAD_LINK="https://zenodo.org/records/4602622/files"

from os.path import join as pjoin

def aggregate_downloaded_batches(wildcards):
    dirsave = checkpoints.download_batches.get(**wildcards).output[0]
    return list(Path(dirsave).rglob("*tar.xz"))

# def aggregate_downloaded_batches_test(wildcards):
#     dirsave = checkpoints.download_batches_test.get(**wildcards).output[0]
#     return list(Path(dirsave).rglob("*tar.xz"))

rule download_verification: 
    input:
        "data/batches_661k_bacteria_all.flag"

rule download_verification_test:
    input:
        "data/batches_661k_bacteria_test.flag"
    
checkpoint download_batches:
    output:
        dirsave=directory("data/bacteria-{subset}")
    input: 
        txt="data/batches_661k_bacteria_{subset}.txt"
    params:
        # dirsave=lambda wildcards: "data+"-{wildcards.subset}",
        link_zenodo=DOWNLOAD_LINK,
    shell:
        """
        mkdir -p {output.dirsave}
        cut {input.txt} -d" " -f3 | head -n 3 | \
        while read f; 
            do wget {params.link_zenodo}/$f -O {output.dirsave}/$f;
        done
        """

rule aggregate_download_batches:
    input:
        aggregate_downloaded_batches
    output:
        flag="data/batches_661k_bacteria_{subset}.flag"
    shell:
        "touch {output.flag}"









# ### --- Test ---

# checkpoint download_batches_test:
#     output:
#         directory(DIRSAVE_TEST)
#     input: 
#         "data/batches_661k_bacteria_test.txt"
#     params:
#         dirsave=DIRSAVE_TEST,
#         link_zenodo=BASEPATH,
#     shell:
#         """
#         mkdir -p {output[0]}
#         cut {input} -d" " -f3 | \
#         while read f; 
#             do wget {params.link_zenodo}/$f -O {params.dirsave}/$f;
#         done
#         """


# rule aggregate_download_batches_test:
#     input:
#         aggregate_downloaded_batches_test
#     output:
#         "data/batches_661k_bacteria_test.flag"
#     shell:
#         "touch {output}"