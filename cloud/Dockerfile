FROM ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends gcc libc6-dev libflint-dev python3 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY dbn.c run.py tasks.json mollifier_y300.txt /app/
RUN gcc -O2 -o /app/dbn /app/dbn.c -lflint -lmpfr -lgmp -lm && /app/dbn block 1000 1100 0.17 0.1414
CMD ["python3", "/app/run.py"]
