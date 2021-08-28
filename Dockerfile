FROM python:3.8.3 

WORKDIR /mirna_profiling_container

COPY . /mirna_profiling_container

RUN pip install -U pip setuptools wheel

RUN pip install --use-feature=in-tree-build --target=/mirna_profiling_container/release .

ENV PYTHONPATH "${PYTHONPATH}:/mirna_profiling_container/release"
ENV PATH "${PATH}:/mirna_profiling_container/release/bin"

ENTRYPOINT ["alignment_stat"]
