git add . 
git commit -m "auto-update"
git push

pip uninstall RiboMetric
pip install git+https://github.com/JackCurragh/RiboMetric.git
