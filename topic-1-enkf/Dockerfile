FROM weiren2/ece408-sp20:amd64
LABEL maintainer "Wei Ren <weiren2@illinois.edu>"

ENV ECE408_BASE /ece408

# tzdata will ask for user input during installation. Disable it.
ARG DEBIAN_FRONTEND=noninteractive

# Topic 1 library dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
		exuberant-ctags \
		libfftw3-dev \
		libgsl-dev \
		liblapacke-dev \
		libopenblas-dev \
		python3 \
		python3-numpy \
		python3-scipy \
		python3-sh \
		python3-matplotlib \
		python3-tk \
		python \
	&& \
	ln -fs /usr/share/zoneinfo/America/Chicago /etc/localtime \
	&& \
	dpkg-reconfigure --frontend noninteractive tzdata \
	&& \
	rm -rf /var/lib/apt/lists/*

### Build Topic 1 files ###

# Make libmdb_matrix
COPY ./libmdb_matrix/ ${ECE408_BASE}/libmdb_matrix
WORKDIR ${ECE408_BASE}/libmdb_matrix
RUN make -j$(nproc) && ls

# Install libconfuse
ADD https://github.com/martinh/libconfuse/releases/download/v3.2/confuse-3.2.tar.gz ${ECE408_BASE}/
WORKDIR ${ECE408_BASE}/
RUN mkdir libconfuse && tar xzf confuse-3.2.tar.gz -C libconfuse --strip-components 1 && rm *.tar.gz
WORKDIR ${ECE408_BASE}/libconfuse 
RUN ls && ./configure && make -j$(nproc) && make install

# Make EnKF - lenkf
COPY ./lenkf/ ${ECE408_BASE}/lenkf
WORKDIR ${ECE408_BASE}/lenkf
ENV PATH_MDB_MATRIX ${ECE408_BASE}/libmdb_matrix
RUN make -j$(nproc)

# Make libmdb_matrix.py visible to lenkf python scripts
RUN ln -s ${ECE408_BASE}/libmdb_matrix/python/libmdb_matrix.py ${ECE408_BASE}/lenkf/python/libmdb_matrix.py
