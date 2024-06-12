#!/bin/bash

# Rstudio server needs write access to /run and /home/rstudio
# Thus, create two writeable directories

mkdir -p $HOME/bind/run $HOME/bind/var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > $HOME/bind/database.conf

FILE=$HOME/GFS/singularity/rocker4csbg.sif

singularity exec \
  --bind $HOME/bind/run:/run \
  --bind $HOME/bind/var-lib-rstudio-server:/var/lib/rstudio-server \
  --bind $HOME/bind/database.conf:/etc/rstudio/database.conf \
  --bind /usr/local/AGFORTELNY:/media/AGFORTELNY \
  --bind /vscratch:/vscratch \
  $FILE rserver \
    --www-address=127.0.0.1 \
    --www-port=8484 \
    --server-user ${USER} \
    --auth-timeout-minutes=0 \
    --auth-stay-signed-in-days=30

