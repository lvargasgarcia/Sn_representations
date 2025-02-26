# Usar la última versión estable de Python como base
FROM ubuntu:latest

# Establecer el directorio de trabajo dentro del contenedor
WORKDIR /Sn_representations

# Copiar la carpeta local al contenedor
COPY . /Sn_representations

# Actualizar repositorios e instalar dependencias del sistema
RUN apt-get update && apt-get install -y \
    gcc \
    libgmp-dev \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/* 

# Comando por defecto al iniciar el contenedor
CMD ["tail", "-f", "/dev/null"]

