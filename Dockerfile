
FROM ubuntu:latest

RUN groupadd -r cimm && useradd -r -g cimm cimm

WORKDIR /home/cimm
ENV FLASK_APP app

RUN apt update && \
    apt install -y build-essential python3.7-dev python3.7-venv python3-venv git && \
    apt clean && \
    python3.7 -m venv venv

LABEL "requirements"="03.06.2019"

COPY requirements.txt boot.sh ./
RUN venv/bin/pip install --no-cache-dir -r requirements.txt && \
    venv/bin/pip install --no-cache-dir gunicorn

COPY DATA DATA
COPY USER USER
WORKDIR /home/cimm/DATA
RUN /home/cimm/venv/bin/pip install -e .
WORKDIR /home/cimm/USER
RUN /home/cimm/venv/bin/pip install -e .

WORKDIR /home/cimm

COPY fragmentor-2017.x bin/
ENV PATH="/home/cimm/bin:${PATH}"

COPY --chown=cimm:cimm zinc.pickle ./
COPY --chown=cimm:cimm app app

USER cimm
EXPOSE 5000
ENTRYPOINT ["./boot.sh"]
