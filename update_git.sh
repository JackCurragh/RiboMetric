git add . 
git commit -m "${1}"
git push

pip uninstall RiboMetric -y
pip install git+https://github.com/JackCurragh/RiboMetric.git

RiboMetric run -b Calviello16/Calviello.transcriptome.sorted.bam -a Calviello16/gencode.v45_RiboMetric.tsv -S 100000 --all --global-offset 15

# RiboMetric run -j test/Calviellotranscriptomesorted_RiboMetric.json --all