FROM openjdk:17-slim-bullseye
#download BBMAP_39.01.tar.gz and place in directory of this Dockerfile
# https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download
#download minimap2-2.24_x64-linux.tar.bz2 and place in directory of this Dockerfile
# https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2

COPY BBMap_39.01.tar.gz ./
COPY minimap2-2.24_x64-linux.tar.bz2 ./

RUN tar -xzvf BBMap_39.01.tar.gz && rm BBMap_39.01.tar.gz

RUN apt-get update && apt-get install -y make bzip2 gcc
RUN cd bbmap/jni && make -f makefile.linux

RUN tar -xjvf minimap2-2.24_x64-linux.tar.bz2 && rm minimap2-2.24_x64-linux.tar.bz2

RUN apt-get install -y python3-pip samtools bcftools

RUN pip3 install --no-cache-dir biopython==1.76 pyfasta pandas pysam
RUN pip3 install --no-cache-dir sklearn numpy matplotlib scipy openpyxl
RUN ln -s /minimap2-2.24_x64-linux/minimap2 /usr/bin/minimap2

ENTRYPOINT ["/bin/bash"]

