FROM quay.io/jgallowa/mypython
RUN python -m pip install --upgrade pip
COPY . .
RUN pip install -r requirements.txt
RUN python setup.py install
