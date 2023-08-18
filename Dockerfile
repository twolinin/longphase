FROM ubuntu:20.04
LABEL maintainer="https://github.com/twolinin/longphase"
LABEL version="1.5.1"

RUN apt-get update && \
    apt-get install -y git g++ gcc autoconf make zlib1g-dev libbz2-dev liblzma-dev && \
    rm -rf /var/lib/apt/lists/* 

WORKDIR /opt/longphase
RUN git clone https://github.com/twolinin/longphase.git /opt/longphase && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase/.git

ENV PATH="${PATH}":${HOME}/bin:/opt/longphase

CMD ["longphase", "phase", "--help"]
#docker run -v "/gpu_disk2/jyunhong104/sup/:/input" -v "/gpu_disk/zhenyu111/longphase/output_longphase/:/output" longphase:1.5.1 longphase phase -s "/input/hg002.wf_snp.vcf.gz" -r "/input/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" -b "/input/hg002.sup.10x.bam" -t 8 -o "/output/longphase" --ont