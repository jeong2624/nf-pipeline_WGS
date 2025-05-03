nextflow run main.nf \
  --input sample_info.csv

rm -rf work .nextflow*
find ./reports -type l -delete
