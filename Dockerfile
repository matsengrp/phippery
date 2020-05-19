FROM python:3.6
RUN python -m pip install --upgrade pip
COPY . .
RUN pip install -r requirements.txt
RUN python setup.py install
