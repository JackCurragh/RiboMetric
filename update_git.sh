git add . 
git commit -m "auto-update"
git push

pip uninstall RiboMetric
pip install git+https://github.com/JackCurragh/RiboMetric.git

RiboMetric run -b Calviello16/Calviello.transcriptome.sorted.bam -a Calviello16/gencode.v45_RiboMetric.tsv -S 1000000 -T 500 --all