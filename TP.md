## TP détection de variants

L'animal SALSA semble présenter un phénotype atypique, une absence de cornes.
Il y a une forte présomption que ce anomalie phénotypique soit d'origine.
Par ailleurs deux autres animaux, ZOUK et TANGO, semble présenter un phénotype intermédiaire, une seul corne.

Pour en avoir le coeur net, on procède à du séquençage WGS de ces animaux pour identifier les variants.
Les différentes étapes ci-dessous vont nous guider dans ce travail.

#### Création de l'environement conda.
```
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/environment.yml
mamba env create -f environment.yml
```

#### Création de l'arborescence
```
mkdir -p data/ref
mkdir -p data/fastq
```
### Récupération ds données
```
cd data/ref
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/reference.fa.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/reference.gff3.bgz
cd ../data/fastq
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/MAMBO_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/MAMBO_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/SALSA_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/SALSA_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/TANGO_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/TANGO_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ZOUK_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ZOUK_R2.fastq.gz
```


