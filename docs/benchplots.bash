#!/bin/bash

# Définir le chemin du dossier benchmark
BENCHMARK_DIR="_data/benchmark"
OUTPUT_DIR="docs/benchplots"

# Vérifier si le dossier benchmark existe
if [ ! -d "$BENCHMARK_DIR" ]; then
  echo "Le dossier $BENCHMARK_DIR n'existe pas."
  exit 1
fi

# Créer le dossier de sortie s'il n'existe pas
mkdir -p "$OUTPUT_DIR"

python ../../tools/workflow_tools/compare_rules.py \
  --benchmarks "$BENCHMARK_DIR" \
  --subsample 1000 \
  --output "$OUTPUT_DIR/compare_rules.png" 

# Parcourir les dossiers dans benchmark
for folder in "$BENCHMARK_DIR"/*; do
  if [ -d "$folder" ]; then
    rule_name=$(basename "$folder")
    
    if [ "$rule_name" == "estimate" ]; then
        extra_params="--size_params index_size_mb"
    else
        extra_params=""
    fi

    echo "Lancement de la commande pour $rule_name..."
    python ../../tools/workflow_tools/plot_benchmark.py \
      --output docs/benchplots/$rule_name.benchmark.png \
      --logx \
      --logy \
      --benchmarks _data/benchmark/$rule_name/*.benchmark.txt \
      --title "RAM and Time performance for rule: $rule_name" $extra_params \
      --regression 0.15 \
      --subsample 500
    
    if [ $? -eq 0 ]; then
      echo "Commande réussie pour $rule_name."
    else
      echo "Erreur lors de l'exécution pour $rule_name."
    fi
  fi
done

echo "Traitement terminé."
