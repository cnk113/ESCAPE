FROM ubuntu:18.04
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' >> /etc/apt/sources.list
RUN apt-get update && apt-get install -y
RUN apt install r-base

RUN pip3 install numpy pandas scipy sklearn tensorflow==1.15rc1 configparser keras torch, networkx, python-louvain

RUN git clone https://github.com/broadinstitute/CellBender.git
RUN pip3 install -e CellBender

RUN git clone https://github.com/AllonKleinLab/scrublet.git
RUN cd scrublet && pip3 install -r requirements.txt && pip3 install --upgrade .

RUN git clone https://github.com/JonathanShor/DoubletDetection.git
RUN cd DoubletDetection && pip3 install .

RUN git clone git://github.com/KrishnaswamyLab/MAGIC.git
RUN cd MAGIC/python && python setup.py install

RUN git clone https://github.com/lanagarmire/deepimpute
RUN cd deepimpute && pip3 install .

RUN wget http://statmath.wu.ac.at/software/RngStreams/rngstreams-1.0.1.tar.gz
RUN tar zxvf rngstreams-1.0.1.tar.gz
RUN cd rngstreams-1.0.1 && ./configure --prefix=/usr/local && make && make install
RUN git clone https://github.com/songfd2018/BUSseq-1.0.git
RUN cd BUSseq-1.0 && make

RUN git clone https://github.com/xikanfeng2/I-Impute.git
RUN cd I-Impute && pip3 install -r requirements.txt && python3 setup.py install

RUN git clone https://github.com/raphael-group/netNMF-sc.git
RUN cd netNMF-sc && python3 setup.py install

RUN git clone https://github.com/AltschulerWu-Lab/scScope.git
RUN cd scScope && python3 setup.py install

RUN pip3 install dca
RUN pip3 install scvi

RUN git clone https://github.com/cnk113/ESCAPE.git
RUN cd ESCAPE Rscript requirements.R
