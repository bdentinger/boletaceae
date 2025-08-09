FROM mambaorg/micromamba:1.5.8
COPY env.yml /tmp/env.yml
RUN micromamba create -y -f /tmp/env.yml -n boletaceae
SHELL ["/bin/bash","-lc"]
WORKDIR /app
COPY . /app
EXPOSE 8080
CMD micromamba run -n boletaceae uvicorn api.main:app --host 0.0.0.0 --port 8080
