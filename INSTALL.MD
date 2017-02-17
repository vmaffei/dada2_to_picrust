# Install Instructions

## Using VirtualEnv
```
virtualenv --no-site-packages env
source env/bin/activate

pip install numpy==1.7.1
pip install scipy
pip install scikit-bio==0.2.3
pip install qiime_default_reference
pip install pynast==1.2.2
cd bin/picrust-1.1.0
python setup.py install
cd ../..
