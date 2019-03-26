#! /bin/bash
# S. Daout 23/02/2012
# Shell facilitant la création du fichier d'entré bridge.in necessaire pour le perl nsb_unwrap.pl

# Initialisations
rm -f bridge.in
rm -f bridge.tmp2
rm -f bridge.tmp

DIRS=./


# Initial checkings
if [ "$1" = "" ]; then
  echo "utilisation : create_bridge.sh bridge.asc"
  echo
  echo "utilisation :"
  echo "mdx.pl filt_date1_date2.int"
  echo "Clic droit sur les ponts"
  echo "Copiez la liste des ponts dans un fichier texte"
  echo "Indiquez la différence de cycle dans la dernière colonne du fichier texte pour chaque pont"
  echo ; exit
fi


fichier=$DIRS/$1

echo 'Extraction des valeurs en azimuth, en range et différence de cycle...'
echo
awk '{print $3,$5,0}' $fichier > bridge.tmp
cat bridge.tmp

echo 'Création du ficher bridge.in'
#gfortran -w -o bridge create_bridge.f95
bridge 
echo 

cat bridge.in
rm bridge.tmp






 

