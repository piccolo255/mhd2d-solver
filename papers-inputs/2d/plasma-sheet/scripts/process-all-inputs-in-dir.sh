#!/bin/bash

set -e

if [[ "$#" = 2 ]]; then
   pushd . >/dev/null 2>&1
   cd "$1"
   INPUTDIR=$PWD
   popd >/dev/null 2>&1
   
   if [[ -d "$2" ]]; then
      read -p "Directory ${2} already exists; continue (y/n)?" -n 1 -r
      echo
      if [[ $REPLY =~ ^[^Yy]$ ]]; then
         echo "Aborted."
         exit 1
      fi
   else
      mkdir -p "$2"
   fi
   
   pushd . >/dev/null 2>&1
   cd "$2"
   DATADIR=$PWD
   popd >/dev/null 2>&1
else
   echo "usage: $0 source_dir data_dir exe_file"
   exit 1
fi

LOGFILE="${DATADIR}/process.log"
SIMEXE=$(readlink -f mhd2d)

echo "$(date +%T) :: Started to process directory '${INPUTDIR}'" | tee $LOGFILE

cd "$INPUTDIR"

for INPUT in *.ini; do
   echo "$(date +%T) :: --- Processing ${INPUT}..." | tee -a $LOGFILE
   
   SIMNAME=$(basename "${INPUT}" .txt)
   SIMDIR="${DATADIR}/${SIMNAME}"
   
   mkdir -p "${SIMDIR}"
   cp "$INPUT" "${SIMDIR}/${INPUT}"
   
   cd "$SIMDIR"
   $SIMEXE "$INPUT" 2>&1 | tee "${INPUT}.log"
   cd "$INPUTDIR"
done

echo "$(date +%T) :: Finished." | tee -a $LOGFILE
