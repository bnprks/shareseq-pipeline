# Setup Notes:

## Building singularity image
- Get Docker Desktop: https://www.docker.com/products/docker-desktop/
- In the repo root directory, run `docker build -t shareseq .`
- After a successful build:
    1. Export docker image to tar file `docker save shareseq > shareseq.tar`
    2. Upload to sherlock
    3. Run `singularity build input.tar output.sif`, recommending to set
       `--disable-cache` and `--tmpdir $L_SCRATCH/singularity`
    - As a one-step script run from the laptop, this is roughly
    ```shell
    export NAME=shareseq && \
    docker save $NAME | ssh nb "
        cat > \$L_SCRATCH/$NAME.tar && 
        mkdir -p \$L_SCRATCH/singularity && 
        bash -lc \"singularity build --disable-cache --force \
            --tmpdir \$L_SCRATCH/singularity \
            \$L_SCRATCH/$NAME.sif \
            docker-archive://\$L_SCRATCH/$NAME.tar\" && 
    cp \$L_SCRATCH/$NAME.sif \$OAK/\$USER/inbox/$NAME.sif
    "
    ```

## Code download links
- Share_seqV2 pipeline from: https://github.com/masai1116/SHARE-seq-alignmentV2/
    - RSeQC download v4.0 tar.gz from: http://rseqc.sourceforge.net/#download-rseqc