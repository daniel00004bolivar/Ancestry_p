mkdir -p out_downsample

for tag in pool1 pool2 hosp; do
  for f in 0.1 0.3 0.5; do
    bam="data/downsample/bam/${tag}_f${f}.sorted.bam"
    out="out_downsample/${tag}_N10_f${f}"
    bash bin/iadmix_master_script.sh \
      --panel_name N10 \
      --chr_list "1-22" \
      --bam "$bam" \
      --output "$out" \
      --make_report
  done
done

mkdir -p out_downsample

for tag in pool1 pool2 hosp; do
  for f in 0.1 0.3 0.5; do
    bam="data/downsample/bam/${tag}_f${f}.sorted.bam"
    out="out_downsample/${tag}_N11_f${f}"
    bash bin/iadmix_master_script.sh \
      --panel_name N11 \
      --chr_list "1-22" \
      --bam "$bam" \
      --output "$out" \
      --make_report
  done
done


